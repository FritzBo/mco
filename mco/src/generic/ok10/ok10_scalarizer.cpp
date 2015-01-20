//
//  ok10_scalarizer.cpp
//  mco
//
//  Created by Fritz Bökler on 16.12.14.
//
//

#include <mco/generic/ok10/ok10_scalarizer.h>

#include <cmath>
#include <iostream>

using std::max;
using std::min;
using std::cout;
using std::endl;

namespace mco {

void Ok10Scalarizer::
Calculate_solutions(std::list<Point *>& solutions) {

    LexStageCompare stage_cmp(epsilon_);
    ComponentwisePointComparator cmp(epsilon_, false);
    EqualityPointComparator eq(epsilon_);

    const Point zero(dimension_);
    Point value(dimension_);
    Point lambda(dimension_);

    std::set<stage_type, LexStageCompare> already_found_stages(stage_cmp);

    dd_set_global_constants();

    pending_stage_queue.clear();


    initial_hyperplane_ = 0;
    double max_ub = -std::numeric_limits<double>::infinity();
    normalization_ = std::numeric_limits<double>::infinity();
    for(double x : upper_bounds_) {
        initial_hyperplane_ += x;
        max_ub = max(max_ub, x);
        normalization_ = min(normalization_, x);
    }
    initial_hyperplane_ *= max_ub;

#ifndef NDEBUG
    cout << "Initial hyperplane offset: " << initial_hyperplane_ << endl;
    cout << "Normalization: " << normalization_ << endl;
#endif

    stage_type initial_stage;
    for(unsigned objective_i = 0; objective_i < dimension_; ++objective_i) {

        Point point(dimension_);
        point[objective_i] = 0x7FFFFFFF;
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

            auto new_point_it = pending_stage.find(value);

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
                    new_stage.insert(value);

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

                auto pos_it = pending_stage.find(test_point);

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
                new_stage.insert(value);

                pending_stage_queue.push_back(new_stage);
            }

        }

    }

#ifndef NDEBUG
    cout << "Total number of stages: " << already_found_stages.size() << endl;
#endif

    dd_free_global_constants();
}

bool Ok10Scalarizer::calculate_lambda(stage_reference stage,
                                      Point& lambda) {

    std::list<Point> hp_defining_vectors;

    dd_MatrixPtr ls = dd_CreateMatrix(dimension_,
                                      dimension_ + 1);

    const Point* ref_point =  &*stage.begin();

    unsigned row_i = 0;

    auto point_it = ++stage.begin();
    while(point_it != stage.end()) {

        for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
            dd_set_si(ls->matrix[row_i][col_i], (int)
                      ref_point->operator[](col_i - 1) - point_it->operator[](col_i - 1));
        }

        dd_set_si(ls->matrix[row_i][0], 0);
        set_addelem(ls->linset, row_i+1);

        ++row_i;
        ++point_it;
    }

    for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
        dd_set_si(ls->matrix[row_i][col_i], 1);
    }
    
    dd_set_si(ls->matrix[row_i][0], -10000);
    set_addelem(ls->linset, row_i+1);
    
    ls->representation = dd_Inequality;
    
#ifndef NDEBUG
    dd_WriteMatrix(stdout, ls);
#endif
    
    dd_ErrorType err;
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(ls, &err);
    
    dd_MatrixPtr sol;
    sol = dd_CopyGenerators(poly);
    
#ifndef NDEBUG
    dd_WriteMatrix(stdout, sol);
#endif
    
    if(sol->rowsize != 1) {
        return false;
    }
    
    for(unsigned i = 1; i < dimension_ + 1; ++i) {
        lambda[i - 1] = dd_get_d(sol->matrix[0][i]);
    }
    
    return true;
}

} // namespace mco
