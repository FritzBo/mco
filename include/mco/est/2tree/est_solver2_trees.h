#pragma once
/*
 * est_solver2_trees.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef EST_SOLVER2_TREES_H_
#define EST_SOLVER2_TREES_H_

#include "../abstract_est_solver.h"

namespace mco {

class ESTSolver2Trees : public mco::AbstractESTSolver {
public:
	ESTSolver2Trees(VectorWeightedGraph &instance) : mco::AbstractESTSolver(instance) {}

	virtual void Solve();
};

} /* namespace mco */
#endif /* EST_SOLVER2_TREES_H_ */
