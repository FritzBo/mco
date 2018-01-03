#pragma once
/*
 * online_vertex_enumerator.h
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#ifndef ONLINE_VERTEX_ENUMERATOR_CDD_H_
#define ONLINE_VERTEX_ENUMERATOR_CDD_H_

#include <vector>
#include <iostream>

#include <cdd/setoper.h>
#include <cdd/cdd_f.h>

#include <mco/basic/point.h>
#include <mco/generic/benson_dual/abstract_online_vertex_enumerator.h>

namespace mco {

class OnlineVertexEnumeratorCDD : AbstractOnlineVertexEnumerator {
public:
	OnlineVertexEnumeratorCDD() = delete;
	virtual ~OnlineVertexEnumeratorCDD();

    bool has_next();
    Point * next_vertex();
    void add_hyperplane(Point& vertex, Point& normal, double rhs);

    unsigned int number_of_hyperplanes() {
		return number_hyperplanes_;
	}

	double get_time() {
//		std::cout << "End cycles: " << cycles_ << std::endl;
//		std::cout << "Time: " << (cycles_/(double) CLOCKS_PER_SEC) << std::endl;
		return cycles_/(double) CLOCKS_PER_SEC;
	}

	OnlineVertexEnumeratorCDD(Point &initial_value, unsigned int dimension, double epsilon);

protected:
	int number_hyperplanes_;

	std::vector<Point> unprocessed_vertices_;

	ddf_MatrixPtr h_representation_;
    unsigned new_line_;

private:
    void expand_h_representation();
};

} /* namespace mco */
#endif /* ONLINE_VERTEX_ENUMERATOR_CDD_H_ */
