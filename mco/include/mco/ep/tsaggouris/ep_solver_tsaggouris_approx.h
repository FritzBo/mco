#pragma once
/*
 * ep_solver_tsaggouris_approx.h
 *
 *  Created on: 26.03.2013
 *      Author: fritz
 */

#ifndef EP_SOLVER_TSAGGOURIS_APPROX_H_
#define EP_SOLVER_TSAGGOURIS_APPROX_H_

#include <vector>

#include <mco/basic/point.h>
#include <mco/ep/basic/abstract_ep_solver.h>
#include <mco/basic/forward_star.h>

namespace mco
{

class EpSolverTsaggourisApprox : public AbstractSolver<std::list<edge>>
{

public:

    struct Instance
    {
        ForwardStar *graph = nullptr;
        FSEdgeArray<Point *> *costs = nullptr;
        unsigned dimension = 0;
        node source = -1;
        node target = -1;
    };

    void Solve(Instance &instance, Point &epsilon);

    ~EpSolverTsaggourisApprox() noexcept override = default;

private:

    struct Label
    {
        const Point point;
        mco::node n;
        const Label *const pred;
        bool mark_dominated;

        inline Label(Point &&point,
                     mco::node n,
                     const Label *pred)
                : point(std::move(point)),
                  n(n),
                  pred(pred),
                  mark_dominated(false)
        {
        }

        inline Label(const Label &label)
                : point(label.point),
                  n(label.n),
                  pred(label.pred),
                  mark_dominated(label.mark_dominated)
        {
        }

        Label &operator=(const Label &label) = delete;
    };

    unsigned position(Instance &instance, const Point &, const Point &) const;

    void ExtendAndMerge(Instance &instance,
                        std::vector<const Label *> &,
                        const mco::edge,
                        const mco::node,
                        std::vector<const Label *> &,
                        const mco::Point &epsilon) const;

    std::vector<double> c_min_;
    std::vector<double> c_max_;
    std::vector<double> max_path_cost_;
    std::vector<double> bases_;

};

} /* namespace mco */
#endif /* EP_SOLVER_TSAGGOURIS_APPROX_H_ */
