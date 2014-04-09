#pragma once
/*
 * abstract_solver.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef ABSTRACT_SOLVER_H_
#define ABSTRACT_SOLVER_H_

#include <list>
#include <iterator>

#include <mco/core/point.h>
#include <mco/core/abstract_graph_instance.h>

namespace mco {

class AbstractSolver {

public:
    
	AbstractSolver() = default;
	AbstractSolver(const AbstractSolver &) = delete;

	AbstractSolver & operator=(const AbstractSolver&) = delete;

	const std::list<const Point *> & solutions() const { return solutions_; }

protected:

	void add_solution(const Point * solution) {
		solutions_.push_back(solution);
	}

	template<class InputIterator>
	void add_solutions(InputIterator begin, InputIterator end) {
		solutions_.insert(solutions_.end(), begin, end);
	}

	void reset_solutions() {
		solutions_.clear();
	}

private:
    
	std::list<const Point *> solutions_;

};
    
}

#endif /* ABSTRACT_SOLVER_H_ */
