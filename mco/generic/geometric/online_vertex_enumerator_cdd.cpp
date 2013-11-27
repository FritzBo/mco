/*
 * online_vertex_enumerator.cpp
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#include <mco/generic/geometric/online_vertex_enumerator_cdd.h>

#include <iostream>
#include <cassert>
#include <list>

using std::list;

namespace mco {

OnlineVertexEnumeratorCDD::~OnlineVertexEnumeratorCDD() {

	dd_FreeMatrix(h_representation_);

	dd_free_global_constants();
}

bool OnlineVertexEnumeratorCDD::has_next() {
	return !unprocessed_vertices_.empty();
}

Point * OnlineVertexEnumeratorCDD::next_vertex() {
	if(unprocessed_vertices_.empty())
		return nullptr;

	Point * p = new Point(unprocessed_vertices_.front());
	unprocessed_vertices_.pop_front();

	return p;
}

void OnlineVertexEnumeratorCDD::add_hyperplane(Point &vertex, Point &normal, double rhs) {
	clock_t start = clock();
//	std::cout << start << std::endl;

//	std::cout << "begin" << std::endl;

//	std::cout << "normal: " << normal << std::endl;
//	std::cout << "rhs: " << rhs << std::endl;

	number_hyperplanes_++;

	dd_PolyhedraPtr poly;
	dd_MatrixPtr v_representation, new_face;
	dd_MatrixPtr new_equality;
	dd_ErrorType err;

	new_equality = dd_CreateMatrix(1, dimension_ + 1);
	for(unsigned int i = 0; i < dimension_; ++i)
		dd_set_d(new_equality->matrix[0][i + 1], normal[i]);
	dd_set_d(new_equality->matrix[0][0], -rhs);

	dd_MatrixAppendTo(&h_representation_, new_equality);

	new_face = dd_CopyMatrix(h_representation_);
	set_addelem(new_face->linset, new_face->rowsize);
	new_face->representation = dd_Inequality;

//	dd_WriteMatrix(stdout, new_face);

	poly = dd_DDMatrix2Poly(new_face, &err);

	if(err != dd_NoError) {
		dd_WriteErrorMessages(stdout, err);
		assert(false);
	}

	v_representation = dd_CopyGenerators(poly);

//	dd_WriteMatrix(stdout, v_representation);

//	std::cout << "removing..." << std::endl;

	auto it = unprocessed_vertices_.begin();
	while(it != unprocessed_vertices_.end()) {
		if(*it * normal - rhs < epsilon_) {
//			std::cout << "Removed " << *it << std::endl;
			it = unprocessed_vertices_.erase(it);
		} else
			it++;
	}

//	std::cout << "adding.." << std::endl;

	for(int i = 0; i < v_representation->rowsize; ++i) {
		if(dd_get_d(v_representation->matrix[i][0]) != 1)
			continue;

		Point p(dimension_);
		for(unsigned j = 0; j < dimension_; ++j)
			p[j] = dd_get_d(v_representation->matrix[i][j + 1]);

		unprocessed_vertices_.push_back(p);

//		std::cout << "Added " << p << std::endl;
	}

	dd_FreeMatrix(v_representation);
	dd_FreeMatrix(new_face);
	dd_FreeMatrix(new_equality);
	dd_FreePolyhedra(poly);

//	std::cout << "end" << std::endl;
//	std::cout << clock() << std::endl;
//	std::cout << clock() - start << std::endl;

	cycles_ += clock() - start;
//	std::cout << cycles_ << std::endl;
}

} /* namespace mco */
