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

class EpSolverWarburtonApprox : public AbstractSolver<std::list<edge>> {

public:
	bool Solve(ogdf::Graph& graph,
               std::function<Point&(edge)> cost_function,
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
};

}

#endif /* EP_SOLVER_WARBURTON_APPROX_H_ */
