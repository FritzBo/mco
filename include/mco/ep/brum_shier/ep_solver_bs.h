#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_H_
#define BSSSA_H_

#include "../abstract_ep_solver.h"

namespace mco {

class EpSolverBS : public mco::AbstractEpSolver {
public:
	explicit EpSolverBS(mco::EpInstance &instance) : AbstractEpSolver(instance) {}
	virtual void Solve();
};

}

#endif /* BSSSA_H_ */
