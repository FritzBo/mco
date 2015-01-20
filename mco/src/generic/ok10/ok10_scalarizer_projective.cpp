//
//  ok10_scalarizer_projective.cpp
//  mco
//
//  Created by Fritz Bökler on 16.12.14.
//
//

#include <mco/generic/ok10/ok10_scalarizer_projective.h>

#include <iostream>
#include <algorithm>

using std::cout;
using std::endl;

namespace mco {

void Ok10ScalarizerProjective::
Calculate_solutions(std::list<Point *>& solutions) {

    LexStageCompare stage_cmp(epsilon_);
    ComponentwisePointComparator cmp(epsilon_, false);
    EqualityPointComparator eq(epsilon_);

    const Point zero(dimension_);
    Point value(dimension_);
    Point projective_value(dimension_ + 1);
    Point lambda(dimension_);

    std::set<stage_type, LexStageCompare> already_found_stages(stage_cmp);

    ddf_set_global_constants();

    pending_stage_queue.clear();

    stage_type initial_stage;
    for(unsigned objective_i = 0; objective_i < dimension_; ++objective_i) {

        Point point(dimension_ + 1);
        point[objective_i] = 1;
        initial_stage.insert(std::move(point));
    }
    pending_stage_queue.push_back(std::move(initial_stage));

    while(!pending_stage_queue.empty()) {

        stage_type pending_stage = pending_stage_queue.front();
        pending_stage_queue.pop_front();

#ifndef NDEBUG
        cout << "Current stage:" << endl;
        for(auto& point : pending_stage) {
            cout << point << endl;
        }
        cout << "Currently " << pending_stage_queue.size() << " stages in queue." << endl;
#endif

        if(already_found_stages.find(pending_stage) != already_found_stages.end()) {
#ifndef NDEBUG
            cout << "Stage already investigated" << endl;
#endif
            continue;
        } else {
            already_found_stages.insert(pending_stage);
        }

        if(!calculate_lambda(pending_stage, lambda)) {
#ifndef NDEBUG
            cout << "Stage is linear dependent." << endl;
#endif
            continue;
        }

#ifndef NDEBUG
        cout << "Lambda: " << lambda << endl;
#endif

        if(cmp(zero, lambda)) {

            solver_(lambda, value);

#ifndef NDEBUG
            cout << "Point found: " << value << endl;
#endif

            for(unsigned i = 0; i < dimension_; ++i) {
                projective_value[i] = value[i];
            }
            projective_value[dimension_] = 1;

            auto new_point_it = pending_stage.find(projective_value);

            // Point is member of the stage
            if(new_point_it != pending_stage.cend()) {

#ifndef NDEBUG
                cout << "Facet found." << endl;
#endif


                // Point is NOT member of the stage
            } else {

                auto pred = [&value, eq] (Point* point) {
                    return eq(&value, point);
                };

                auto it = std::find_if(solutions.begin(),
                                       solutions.end(),
                                       pred);

                if(it == solutions.end()) {
                    solutions.push_back(new Point(value));
#ifndef NDEBUG
                    cout << "Point is new extreme point." << endl;
#endif

                }

#ifndef NDEBUG
                cout << "Creating new stages." << endl;
#endif

                for(unsigned i = 0; i < dimension_; ++i) {
                    stage_type new_stage;

                    unsigned j = 0;
                    auto old_point_it = pending_stage.begin();
                    while(old_point_it != pending_stage.end()) {
                        if(i != j) {
                            new_stage.insert(*old_point_it);
                        }

                        ++old_point_it;
                        ++j;
                    }
                    new_stage.insert(projective_value);

                    pending_stage_queue.push_back(new_stage);
                }

            }

        } else {

#ifndef NDEBUG
            cout << "Negative components. Creating help stages." << endl;
#endif

            auto ex_point_it = solutions.begin();
            while(ex_point_it != solutions.end()) {
                Point& test_point = **ex_point_it;

                for(unsigned i = 0; i < dimension_; ++i) {
                    projective_value[i] = test_point[i];
                }
                projective_value[dimension_] = 1;

                auto pos_it = pending_stage.find(projective_value);

                if(pos_it == pending_stage.end()) {
                    value = test_point;
                    break;
                }

                ++ex_point_it;

            }

            if(ex_point_it == solutions.end()) {
                pending_stage_queue.push_back(pending_stage);
                continue;
            }

#ifndef NDEBUG
            cout << "Dummy Point: " << **ex_point_it << endl;
#endif

            for(unsigned i = 0; i < dimension_; ++i) {
                stage_type new_stage;

                unsigned j = 0;
                auto old_point_it = pending_stage.begin();
                while(old_point_it != pending_stage.end()) {
                    if(i != j) {
                        new_stage.insert(*old_point_it);
                    }

                    ++old_point_it;
                    ++j;
                }
                new_stage.insert(projective_value);

                pending_stage_queue.push_back(new_stage);
            }

        }

    }

#ifndef NDEBUG
    cout << "Total number of stages: " << already_found_stages.size() << endl;
#endif

    ddf_free_global_constants();
}

bool Ok10ScalarizerProjective::calculate_lambda(stage_reference stage,
                                                Point& lambda) {

    EqualityPointComparator eq(1E-9);

    std::list<Point> hp_defining_vectors;

    const Point* ref_point = nullptr;
    auto point_it = stage.begin();
    while(point_it != stage.end()) {
        if(point_it->operator[](dimension_) != 0) {
            ref_point = &*point_it;
            break;
        }
        ++point_it;
    }

    if(ref_point == nullptr) {
        ref_point = &*stage.begin();
    }

    unsigned row_i = 0;

    ddf_MatrixPtr ls = ddf_CreateMatrix(dimension_,
                                        dimension_ + 1);

    point_it = stage.begin();
    while(point_it != stage.end()) {

        if(!eq(*ref_point, *point_it)) {

            // Point is an extreme point or all points are rays
            if(ref_point->operator[](dimension_) == 0 || point_it->operator[](dimension_) != 0) {

                for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
                    ddf_set_d(ls->matrix[row_i][col_i], (int)
                              ref_point->operator[](col_i - 1) - point_it->operator[](col_i - 1));
                }

                ddf_set_d(ls->matrix[row_i][0], 0);
                set_addelem(ls->linset, row_i + 1);

                // Point represents an extreme ray and there is at least
                // one extreme point
            } else {

                for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
                    ddf_set_d(ls->matrix[row_i][col_i], point_it->operator[](col_i - 1));
                }

                ddf_set_d(ls->matrix[row_i][0], 0);
                set_addelem(ls->linset, row_i + 1);
            }
            
            ++row_i;
        }
        
        ++point_it;
    }
    
    for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
        ddf_set_si(ls->matrix[row_i][col_i], 1);
    }
    
    ddf_set_si(ls->matrix[row_i][0], -1);
    set_addelem(ls->linset, row_i+1);
    
    ls->representation = ddf_Inequality;
    
#ifndef NDEBUG
    ddf_WriteMatrix(stdout, ls);
#endif
    
    ddf_ErrorType err;
    ddf_PolyhedraPtr poly = ddf_DDMatrix2Poly(ls, &err);
    
    ddf_MatrixPtr sol;
    sol = ddf_CopyGenerators(poly);
    
#ifndef NDEBUG
    ddf_WriteMatrix(stdout, sol);
#endif
    
    if(sol->rowsize != 1) {
        return false;
    }
    
    for(unsigned i = 1; i < dimension_ + 1; ++i) {
        lambda[i - 1] = ddf_get_d(sol->matrix[0][i]);
    }
    
    return true;
}

}
