#pragma once
/*
 * dual_benson_scalarizer.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef DUAL_BENSON_SCALARIZER_H_
#define DUAL_BENSON_SCALARIZER_H_

#include <list>

#include <mco/point.h>
#include <mco/generic/benson_weightspace/dual_benson_vertex_container.h>

namespace mco {

class ScalarizationSolver {
public:
	virtual double Solve_scalarization(Point &weights, Point &value) = 0;
	virtual ~ScalarizationSolver() {
	}
};

class DualBensonScalarizer {
public:
	DualBensonScalarizer(ScalarizationSolver &solver, unsigned int dimension, double epsilon = 1E-6);
	virtual ~DualBensonScalarizer() {
	}

	void Calculate_solutions(std::list<Point *> &solutions);

	double vertex_enumeration_time() {
		if(vertex_container == nullptr)
			return 0;

		return vertex_container->get_time();
	}

	int number_vertices() {
		return vertices_;
	}

	int number_facets() {
		return facets_;
	}

protected:
	unsigned int dimension_;
	double epsilon_;

	ScalarizationSolver &solver_;

private:
	DualBensonVertexContainer *vertex_container;

	int vertices_;
	int facets_;
};

} /* namespace mco */
#endif /* DUAL_BENSON_SCALARIZER_H_ */
