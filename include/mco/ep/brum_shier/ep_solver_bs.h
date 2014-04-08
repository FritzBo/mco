#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_H_
#define BSSSA_H_

#include "../basic/abstract_ep_solver.h"

namespace mco {

class EpSolverBS : public AbstractEpSolver {
public:
	explicit EpSolverBS(EpInstance &instance, double epsilon) : AbstractEpSolver(instance), epsilon_(epsilon) {}
	virtual void Solve();
    
private:
    const double epsilon_;
};

}

#endif /* BSSSA_H_ */
