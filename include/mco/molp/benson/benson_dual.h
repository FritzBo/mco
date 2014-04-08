#pragma once
/*
 * benson_dual.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef BENSON_DUAL_H_
#define BENSON_DUAL_H_

#include <mco/core/abstract_solver.h>

namespace mco {

class DualBensonMolpSolver : AbstractSolver {
public:
	DualBensonMolpSolver();

	void Solve();
};

} /* namespace mco */
#endif /* BENSON_DUAL_H_ */
