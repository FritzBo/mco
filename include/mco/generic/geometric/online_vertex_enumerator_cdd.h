#pragma once
/*
 * online_vertex_enumerator.h
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#ifndef ONLINE_VERTEX_ENUMERATOR_CDD_H_
#define ONLINE_VERTEX_ENUMERATOR_CDD_H_

#include <list>
#include <iostream>

#include <setoper.h>
#include <cdd.h>

#include <mco/point.h>

namespace mco {

class OnlineVertexEnumeratorCDD {
public:
	OnlineVertexEnumeratorCDD() = delete;
	virtual ~OnlineVertexEnumeratorCDD();

	virtual bool has_next();
	virtual Point * next_vertex();
	virtual void add_hyperplane(Point& vertex, Point& normal, double rhs);

	virtual unsigned int number_of_hyperplanes() {
		dd_ErrorType err;
		dd_rowset redundant_rows = dd_RedundantRows(h_representation_, &err);
		return h_representation_->rowsize - set_card(redundant_rows);
	}

	double get_time() {
//		std::cout << "End cycles: " << cycles_ << std::endl;
//		std::cout << "Time: " << (cycles_/(double) CLOCKS_PER_SEC) << std::endl;
		return cycles_/(double) CLOCKS_PER_SEC;
	}


protected:
	OnlineVertexEnumeratorCDD(unsigned int dimension, double epsilon) :
	dimension_(dimension), epsilon_(epsilon), number_hyperplanes_(0), cycles_(0), h_representation_(nullptr) {
		dd_set_global_constants();
	}

	const unsigned int dimension_;
	const double epsilon_;

	int number_hyperplanes_;

	std::list<Point> unprocessed_vertices_;

	clock_t cycles_;

	dd_MatrixPtr h_representation_;
};

} /* namespace mco */
#endif /* ONLINE_VERTEX_ENUMERATOR_CDD_H_ */
