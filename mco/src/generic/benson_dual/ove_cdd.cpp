/*
 * online_vertex_enumerator.cpp
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#include <mco/generic/benson_dual/ove_cdd.h>

#include <iostream>
#include <cassert>
#include <list>

using std::list;

namespace mco {

OnlineVertexEnumeratorCDD::OnlineVertexEnumeratorCDD(Point &initial_value, unsigned int dimension, double epsilon) :
	AbstractOnlineVertexEnumerator(dimension, epsilon),
	number_hyperplanes_(0),
	h_representation_(nullptr) {

	ddf_set_global_constants();

	h_representation_ = ddf_CreateMatrix(2 * (dimension + 1), dimension + 1);

	for(unsigned int i = 0; i < dimension_ - 1; ++i) {
		for(unsigned int j = 0; j < dimension_; ++j)
			ddf_set_d(h_representation_->matrix[i][j + 1], i == j ? 1 : 0);
		ddf_set_d(h_representation_->matrix[i][0], 0);

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
		ddf_set_d(h_representation_->matrix[dimension_ - 1][j + 1], -1);
	ddf_set_d(h_representation_->matrix[dimension_ - 1][dimension_], 0);
	ddf_set_d(h_representation_->matrix[dimension_ - 1][0], 1);

	for(unsigned int j = 0; j < dimension_ - 1; ++j)
		ddf_set_d(h_representation_->matrix[dimension_][j + 1], initial_value[j] - initial_value[dimension_ - 1]);
	ddf_set_d(h_representation_->matrix[dimension_][dimension_], -1);
	ddf_set_d(h_representation_->matrix[dimension_][0], initial_value[dimension_ - 1]);

	h_representation_->representation = ddf_Inequality;
    new_line_ = dimension + 1;

//	dd_WriteMatrix(stdout, h_representation_);
//
//	for(Point &p : unprocessed_vertices_)
//		cout << p << endl;
}

OnlineVertexEnumeratorCDD::~OnlineVertexEnumeratorCDD() {

	ddf_FreeMatrix(h_representation_);

	ddf_free_global_constants();
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

	ddf_PolyhedraPtr poly;
	ddf_MatrixPtr v_representation;
	ddf_ErrorType err;

//    if(new_line_ == dimension_ + 1) {
//        cout << ddf_almostzero << endl;
//    }

    if(new_line_ >= h_representation_->rowsize) {
        expand_h_representation();
    }

	for(unsigned int i = 0; i < dimension_; ++i) {
		ddf_set_d(h_representation_->matrix[new_line_][i + 1], normal[i]);
    }
	ddf_set_d(h_representation_->matrix[new_line_][0], -rhs);

//	new_face = ddf_CopyMatrix(h_representation_);
//	set_addelem(new_face->linset, new_face->rowsize);
//	new_face->representation = ddf_Inequality;

    set_addelem(h_representation_->linset, new_line_);

    h_representation_->representation = ddf_Inequality;
    h_representation_->numbtype = ddf_Real;

//	ddf_WriteMatrix(stdout, h_representation_);

	poly = ddf_DDMatrix2Poly2(h_representation_, ddf_MaxCutoff,&err);

	if(err != ddf_NoError) {
		ddf_WriteErrorMessages(stdout, err);
		assert(false);
        exit(-1);
	}

	v_representation = ddf_CopyGenerators(poly);

    set_delelem(h_representation_->linset, new_line_);

//	ddf_WriteMatrix(stdout, h_representation_);
//    ddf_WriteMatrix(stdout, v_representation);

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
		if(ddf_get_d(v_representation->matrix[i][0]) != 1)
			continue;

		Point p(dimension_);
		for(unsigned j = 0; j < dimension_; ++j)
			p[j] = ddf_get_d(v_representation->matrix[i][j + 1]);

		unprocessed_vertices_.push_back(p);

//		std::cout << "Added " << p << std::endl;
	}

    ++new_line_;

	ddf_FreeMatrix(v_representation);
	ddf_FreePolyhedra(poly);

//	std::cout << "end" << std::endl;
//	std::cout << clock() << std::endl;
//	std::cout << clock() - start << std::endl;

	cycles_ += clock() - start;
//	std::cout << cycles_ << std::endl;
}

void OnlineVertexEnumeratorCDD::expand_h_representation() {
    ddf_MatrixPtr new_h_representation = ddf_CreateMatrix(2 * new_line_, dimension_ + 1);

    for(unsigned row = 0; row < new_line_; ++row) {
        for(unsigned col = 0; col < dimension_ + 1; ++col) {
            ddf_set_d(new_h_representation->matrix[row][col],
                      ddf_get_d(h_representation_->matrix[row][col]));
        }
    }

    ddf_FreeMatrix(h_representation_);
    h_representation_ = new_h_representation;
}

} /* namespace mco */
