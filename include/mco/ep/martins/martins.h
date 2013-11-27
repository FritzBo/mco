#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef MARTINS_B_H_
#define MARTINS_B_H_

#include "../abstract_ep_solver.h"

namespace mco {

class EpSolverMartins : public AbstractEpSolver {

public:
	explicit EpSolverMartins(mco::EpInstance &instance);
	virtual void Solve();

};

}

#endif /* MARTINS_H_ */
