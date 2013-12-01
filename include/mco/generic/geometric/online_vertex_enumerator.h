#pragma once
/*
 * online_vertex_enumerator.h
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#ifndef ONLINE_VERTEX_ENUMERATOR_H_
#define ONLINE_VERTEX_ENUMERATOR_H_

#include <vector>
#include <map>
#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/point.h>

namespace mco {

class OnlineVertexEnumerator {
public:
	OnlineVertexEnumerator() = delete;
	virtual ~OnlineVertexEnumerator();

	virtual bool has_next();
	virtual Point * next_vertex();
	virtual void add_hyperplane(Point &vertex, Point &normal, double rhs);

	virtual unsigned int number_of_hyperplanes();

	double get_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}


protected:
	OnlineVertexEnumerator(unsigned int dimension, double epsilon) :
		dimension_(dimension),
		epsilon_(epsilon),
		node_points_(vertex_graph_, nullptr),
		comp_(epsilon),
		point_nodes_(comp_),
		node_inequality_indices_(vertex_graph_, nullptr),
		cycles_(0) {
	}

	const unsigned int dimension_;
	const double epsilon_;
	ogdf::Graph vertex_graph_;
	ogdf::NodeArray<Point *> node_points_;

	LexicographicPointComparator comp_;
	std::map<Point *, ogdf::node, LexicographicPointComparator> point_nodes_;
	std::list<Point *> unprocessed_projective_points_;

	std::list<Point *> list_of_inequalities_;
	ogdf::NodeArray<std::list<int> *> node_inequality_indices_;

	ogdf::node get_node(Point &non_projective_point);
	Point * to_projective(Point &non_projective_point);
	Point normalize_projective(Point projective_point);

	bool inside_face(Point &p1, Point &p2, bool nondegenerate);

	clock_t cycles_;
};

} /* namespace mco */
#endif /* ONLINE_VERTEX_ENUMERATOR_H_ */
