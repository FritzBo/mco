#pragma once
/*
 * abstract_est_solver.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_EST_SOLVER_H_
#define ABSTRACT_EST_SOLVER_H_

#include <list>

#include "../abstract_solver.h"
#include "../vector_weighted_graph.h"

namespace mco {

class AbstractESTSolver : public AbstractMOSolver {
	VectorWeightedGraph & instance_;
	std::list<Point *> solutions_;

protected:

	VectorWeightedGraph & instance() {
		return instance_;
	}

public:
	AbstractESTSolver() = delete;
	AbstractESTSolver(VectorWeightedGraph & instance) : instance_(instance) {}

};

} /* namespace mco */
#endif /* ABSTRACT_EST_SOLVER_H_ */
