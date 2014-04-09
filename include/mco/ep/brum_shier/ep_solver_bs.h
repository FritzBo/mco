#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_H_
#define BSSSA_H_

#include <mco/core/abstract_solver.h>

namespace mco {

class EpSolverBS : public AbstractSolver {
    
public:
	EpSolverBS(double epsilon = 1E-8)
    : epsilon_(epsilon) { }
    
	virtual void Solve(ogdf::Graph graph,
                       std::function<const Point*(const ogdf::edge)> costs,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target);
    
private:
    const double epsilon_;
};

}

#endif /* BSSSA_H_ */
