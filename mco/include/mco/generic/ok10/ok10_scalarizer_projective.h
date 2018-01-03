//
//  ok10_scalarizer.h
//  mco
//
//  Created by Fritz Bökler on 16.12.14.
//
//

#ifndef mco_ok10_scalarizer_projective_h
#define mco_ok10_scalarizer_projective_h

#include <list>
#include <functional>
#include <set>

#include <cdd/setoper.h>
#include <cdd/cdd.h>

#include <mco/basic/point.h>

namespace mco {


class Ok10ScalarizerProjective {
public:
    Ok10ScalarizerProjective(std::function<double(const Point& weighting, Point& value)> solver,
                   Point& upper_bounds,
                   unsigned int dimension)
    :   dimension_(dimension),
        epsilon_(1E-8),
        upper_bounds_(upper_bounds),
        solver_(solver) {

        assert(dimension_ > 1);
    }

    void Calculate_solutions(std::list<Point *>& solutions);

protected:
    unsigned int dimension_;
    double epsilon_;
    Point& upper_bounds_;

    double initial_hyperplane_;
    double normalization_;

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


} /* namespace mco */

#endif /* namespace mco_ok10_scalarizer_projective_h */
