//
//  ok10_scalarizer.h
//  mco
//
//  Created by Fritz Bökler on 14.11.14.
//
//

#ifndef mco_ok10_scalarizer_h
#define mco_ok10_scalarizer_h

#include <list>
#include <functional>

#include <setoper.h>
#include <cdd.h>

#include <mco/basic/point.h>

namespace mco {


class Ok10Scalarizer {
public:
    Ok10Scalarizer(std::function<double(const Point& weighting, Point& value)> solver,
                   unsigned int dimension)
    :   dimension_(dimension),
    epsilon_(1E-8),
    solver_(solver) {
        assert(dimension_ > 1);
    }

    void Calculate_solutions(std::list<Point *>& solutions);

protected:
    unsigned int dimension_;
    double epsilon_;

    std::function<double(const Point&, Point&)> solver_;

private:
    using stage_type = std::set<Point, LexPointComparator>;
    using stage_pointer = stage_type *;
    using stage_reference = stage_type &;

    std::list<stage_type> pending_stage_queue;

    bool calculate_lambda(stage_reference stage,
                          Point& lambda);

    class LexStageCompare {
    public:
        LexStageCompare(double epsilon = 1E-8)
        :   lex_cmp(epsilon) {}

        bool operator()(const stage_type& stage1,
                        const stage_type& stage2) {

            auto it1 = stage1.begin();
            auto it2 = stage2.begin();
            while(it1 != stage1.end()) {
                if(lex_cmp(*it1, *it2)) {
                    return true;
                } else if(lex_cmp(*it2, *it1)) {
                    return false;
                }
                ++it1; ++it2;
            }
            return false;
        }
    private:
        LexPointComparator lex_cmp;
    };

};

void Ok10Scalarizer::
Calculate_solutions(std::list<Point *>& solutions) {

    LexStageCompare stage_cmp(epsilon_);
    ComponentwisePointComparator cmp(epsilon_);
    EqualityPointComparator eq(epsilon_);

    const Point zero(dimension_);
    Point value(dimension_);
    Point lambda(dimension_);

    std::set<stage_type, LexStageCompare> already_found_stages(stage_cmp);

    ddf_set_global_constants();

    pending_stage_queue.clear();

    stage_type initial_stage;
    for(unsigned objective_i = 0; objective_i < dimension_; ++objective_i) {

        Point point(dimension_);
        point[objective_i] = 1E5;
        initial_stage.insert(std::move(point));
    }
    pending_stage_queue.push_back(std::move(initial_stage));

    while(!pending_stage_queue.empty()) {

        stage_type pending_stage = pending_stage_queue.back();
        pending_stage_queue.pop_back();

        if(already_found_stages.find(pending_stage) != already_found_stages.end()) {
#ifndef NDEBUG
            cout << "Stage already investigated" << endl;
#endif
            continue;
        } else {
            already_found_stages.insert(pending_stage);
        }

#ifndef NDEBUG
        cout << "Current stage:" << endl;
        for(auto& point : pending_stage) {
            cout << point << endl;
        }
        cout << "Currently " << pending_stage_queue.size() << " stages in queue." << endl;
#endif


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

            Point old_ex(dimension_);

            auto ex_point_it = solutions.begin();
            while(ex_point_it != solutions.end()) {
                Point& test_point = **ex_point_it;

                auto pos_it = pending_stage.find(test_point);

                if(pos_it == pending_stage.end()) {
                    old_ex = test_point;
                    break;
                }

                ++ex_point_it;

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
                new_stage.insert(old_ex);

                pending_stage_queue.push_back(new_stage);
            }

        }

    }

//#ifndef NDEBUG
    cout << "Total number of stages: " << already_found_stages.size() << endl;
//#endif

    ddf_free_global_constants();
}

bool Ok10Scalarizer::calculate_lambda(stage_reference stage,
                                      Point& lambda) {

    std::list<Point> hp_defining_vectors;

    ddf_MatrixPtr ls = ddf_CreateMatrix(dimension_,
                                        dimension_ + 1);

    const Point* ref_point =  &*stage.begin();

    unsigned row_i = 0;

    auto point_it = ++stage.begin();
    while(point_it != stage.end()) {

        for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
            ddf_set_si(ls->matrix[row_i][col_i], (int)
                      ref_point->operator[](col_i - 1) - point_it->operator[](col_i - 1));
        }

        ddf_set_si(ls->matrix[row_i][0], 0);
        set_addelem(ls->linset, row_i+1);

        ++row_i;
        ++point_it;
    }

    for(unsigned col_i = 1; col_i < dimension_ + 1; ++col_i) {
        ddf_set_si(ls->matrix[row_i][col_i], 1);
    }

    ddf_set_si(ls->matrix[row_i][0], -100);
    set_addelem(ls->linset, row_i+1);

    ls->representation = ddf_Inequality;

//    ddf_WriteMatrix(stdout, ls);

    ddf_ErrorType err;
    ddf_PolyhedraPtr poly = ddf_DDMatrix2Poly(ls, &err);

    ddf_MatrixPtr sol;
    sol = ddf_CopyGenerators(poly);

//    ddf_WriteMatrix(stdout, sol);

    if(sol->rowsize != 1) {
        return false;
    }

    for(unsigned i = 1; i < dimension_ + 1; ++i) {
        lambda[i - 1] = ddf_get_d(sol->matrix[0][i]);
    }

    return true;
}

} /* namespace mco */

#endif
