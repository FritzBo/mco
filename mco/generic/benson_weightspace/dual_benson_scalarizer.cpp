/*
 * dual_benson_scalarizer.cpp
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#include <mco/generic/benson_weightspace/dual_benson_scalarizer.h>

#include <iostream>
#include <list>

using std::list;

namespace mco {

DualBensonScalarizer::DualBensonScalarizer(ScalarizationSolver &solver, unsigned int dimension, double epsilon) :
	dimension_(dimension), epsilon_(epsilon), solver_(solver), vertex_container(nullptr), vertices_(0), facets_(0) {
}

void DualBensonScalarizer::Calculate_solutions(list<Point *> &solutions) {
	if(this == nullptr) return;

	int nondominated_values = 1;
	int iteration_counter = 0;
	int weighting_counter = 1;

	Point v(dimension_);
	Point value(dimension_);

	for(unsigned int i = 0; i < dimension_ - 1; ++i)
		v[i] = 0;
	v[0] = 1;

	v[dimension_ - 1] = solver_.Solve_scalarization(v, value);

	vertex_container = new DualBensonVertexContainer(value, dimension_, epsilon_);
	vertex_container->next_vertex();

	Point *candidate, weighting(dimension_), inequality(dimension_);
	double scalar_value;
	while(vertex_container->has_next()) {
		iteration_counter++;

//		std::cout << "Iteration: " << iteration_counter << std::endl;
		candidate = vertex_container->next_vertex();

//		std::cout << "New candidate: " << *candidate << std::endl;

		double sum = 0;
		for(unsigned int i = 0; i < dimension_ - 1; ++i) {
			weighting[i] = (*candidate)[i];
			sum += (*candidate)[i];
		}
		weighting[dimension_ - 1] = 1 - sum;

		for(unsigned int i = 0; i < dimension_; ++i)
			value[i] = 0;

//		std::cout << "weighting: " << weighting << std::endl;

		scalar_value = solver_.Solve_scalarization(weighting, value);

//		std::cout << "scalar value: " << scalar_value << std::endl;
//		std::cout << "value vector: " << value << std::endl;

		if(scalar_value - (*candidate)[dimension_ - 1] > -epsilon_) {
			weighting_counter++;
			continue;
		}

		for(unsigned int i = 0; i < dimension_ - 1; ++i)
			inequality[i] = value[i] - value[dimension_ - 1];
		inequality[dimension_ - 1] = -1;

		vertex_container->add_hyperplane(*candidate, inequality, -value[dimension_ - 1]);
		nondominated_values++;

		solutions.push_back(new Point(value));

	}

	vertices_ = nondominated_values;
	facets_ = weighting_counter;

//	std::cout << "Found " << nondominated_values << " nondominated value vectors in " << iteration_counter << " iterations." << std::endl;
//	std::cout << "Where " << weighting_counter << " weightings have been explored." << std::endl;

	delete vertex_container;
}

} /* namespace mco */
