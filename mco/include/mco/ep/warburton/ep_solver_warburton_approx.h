#pragma once
/*
 * ep_solver_warburton_approx.h
 *
 *  Created on: 02.04.2013
 *      Author: fritz
 */

#ifndef EP_SOLVER_WARBURTON_APPROX_H_
#define EP_SOLVER_WARBURTON_APPROX_H_

#include <mco/ep/basic/abstract_ep_solver.h>
#include <mco/ep/basic/ep_instance.h>

namespace mco {

class EpSolverWarburtonApprox : public AbstractSolver<std::list<ogdf::edge>> {

public:
	bool Solve(ogdf::Graph& graph,
               std::function<Point&(ogdf::edge)> cost_function,
               unsigned dimension,
               const ogdf::node source,
               const ogdf::node target,
               const Point& epsilon,
               const Point& bound,
               bool test_only = false,
               bool directed = false,
               unsigned int processes = 2,
               double theta = 2.0);

	~EpSolverWarburtonApprox() noexcept {}
    
private:
    bool computeBounds(const ogdf::Graph& graph,
                       std::function<Point&(ogdf::edge)> cost_function,
                       const Point& epsilon,
                       const Point& bound,
                       const ogdf::node source,
                       const ogdf::node target,
                       double theta,
                       unsigned dimension,
                       bool directed,
                       Point& lb,
                       Point& ub,
                       Point& label_limits);
    
    unsigned compute_skip_function(unsigned dimension,
                                   const Point& lb,
                                   const Point& ub);
    
    bool lagrange_prune(const Graph& graph,
                        std::function<Point*(ogdf::edge)> cost_function,
                        const ogdf::node source,
                        const ogdf::node target,
                        bool directed,
                        Point label_limits,
                        unsigned skip_function,
                        unsigned dimension);
};

}

#endif /* EP_SOLVER_WARBURTON_APPROX_H_ */
