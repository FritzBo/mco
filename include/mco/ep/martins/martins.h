#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef MARTINS_B_H_
#define MARTINS_B_H_

#include <mco/ep/basic/abstract_ep_solver.h>

namespace mco {

class EpSolverMartins : public AbstractEpSolver {

public:
	explicit EpSolverMartins(mco::EpInstance &instance,
                             double epsilon = 0)
    : AbstractEpSolver(instance), epsilon_(epsilon) { }
    
	virtual void Solve();
    
private:
    const double epsilon_;

};

}

#endif /* MARTINS_H_ */
