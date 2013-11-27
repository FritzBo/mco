#pragma once
/*
 * AbstractAPSolver.h
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#ifndef ABSTRACTAPSOLVER_H_
#define ABSTRACTAPSOLVER_H_

#include <memory>

#include <mco/abstract_solver.h>
#include <mco/ap/assignment_instance.h>

namespace mco {

class AbstractAPSolver : public AbstractMOSolver {
public:
	AbstractAPSolver(std::shared_ptr<AssignmentInstance> instance) :
	instance_(instance) {}

protected:
	AssignmentInstance &instance() {
		return *instance_;
	}

private:
	std::shared_ptr<AssignmentInstance> instance_;

};
}


#endif /* ABSTRACTAPSOLVER_H_ */
