/*
 * dual_benson_vertex_container.cpp
 *
 *  Created on: 01.10.2013
 *      Author: fritz
 */

#include <mco/generic/benson_weightspace/dual_benson_vertex_container_cdd.h>

using std::make_pair;

#include <ogdf/basic/Graph.h>

using ogdf::node;

#include <mco/point.h>

namespace mco {

DualBensonVertexContainerCDD::DualBensonVertexContainerCDD(Point &initial_value, unsigned int dimension, double epsilon) :
	OnlineVertexEnumeratorCDD(dimension, epsilon) {

	h_representation_ = dd_CreateMatrix(dimension + 1, dimension + 1);

	for(unsigned int i = 0; i < dimension_ - 1; ++i) {
		for(unsigned int j = 0; j < dimension_; ++j)
			dd_set_d(h_representation_->matrix[i][j + 1], i == j ? 1 : 0);
		dd_set_d(h_representation_->matrix[i][0], 0);

		Point p(dimension_);
		for(unsigned int j = 0; j < dimension_ - 1; ++j)
			p[j] = i == j ? 1 : 0;
		p[dimension_ - 1] = initial_value[i];

		unprocessed_vertices_.push_back(p);
	}

	Point p(dimension);
	p[dimension_ - 1] = initial_value[dimension_ - 1];

	unprocessed_vertices_.push_back(p);

	for(unsigned int j = 0; j < dimension_ - 1; ++j)
		dd_set_d(h_representation_->matrix[dimension_ - 1][j + 1], -1);
	dd_set_d(h_representation_->matrix[dimension_ - 1][dimension_], 0);
	dd_set_d(h_representation_->matrix[dimension_ - 1][0], 1);

	for(unsigned int j = 0; j < dimension_ - 1; ++j)
		dd_set_d(h_representation_->matrix[dimension_][j + 1], initial_value[j] - initial_value[dimension_ - 1]);
	dd_set_d(h_representation_->matrix[dimension_][dimension_], -1);
	dd_set_d(h_representation_->matrix[dimension_][0], initial_value[dimension_ - 1]);

	h_representation_->representation = dd_Inequality;

//	dd_WriteMatrix(stdout, h_representation_);
//
//	for(Point &p : unprocessed_vertices_)
//		cout << p << endl;
}

} /* namespace mco */
