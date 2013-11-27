#pragma once
/*
 * dual_benson_scalarizer.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef DUAL_BENSON_SCALARIZER_CDD_H_
#define DUAL_BENSON_SCALARIZER_CDD_H_

#include <list>

#include <mco/point.h>
#include <mco/generic/benson_weightspace/dual_benson_vertex_container_cdd.h>

namespace mco {

class ScalarizationSolverCDD {
public:
	virtual double Solve_scalarization(Point &weights, Point &value) = 0;
	virtual ~ScalarizationSolverCDD() {
	}
};

class DualBensonScalarizerCDD {
public:
	DualBensonScalarizerCDD(ScalarizationSolverCDD &solver, unsigned int dimension, double epsilon = 1E-6);
	virtual ~DualBensonScalarizerCDD() {
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

	ScalarizationSolverCDD &solver_;

private:
	DualBensonVertexContainerCDD *vertex_container;

	int vertices_;
	int facets_;
};

} /* namespace mco */
#endif /* DUAL_BENSON_SCALARIZER_CDD_H_ */
