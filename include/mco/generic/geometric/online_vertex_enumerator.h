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
#include <queue>
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
		unprocessed_projective_points_(LexicographicPointComparator(epsilon, false)),
		node_inequality_indices_(vertex_graph_, nullptr),
		birth_index_(vertex_graph_, -1),
		cycles_(0) {
	}

	const unsigned int dimension_;
	const double epsilon_;
	ogdf::Graph vertex_graph_;
	ogdf::NodeArray<Point *> node_points_;

	LexicographicPointComparator comp_;
	std::map<Point *, ogdf::node, LexicographicPointComparator> point_nodes_;
	std::priority_queue<Point *, std::vector<Point *>, LexicographicPointComparator> unprocessed_projective_points_;

	std::vector<Point *> list_of_inequalities_;
	ogdf::NodeArray<std::list<int> *> node_inequality_indices_;
	ogdf::NodeArray<unsigned int> birth_index_;

	ogdf::node get_node(Point &non_projective_point);
	Point * to_projective(Point &non_projective_point);
	Point normalize_projective(Point projective_point);

	bool inside_face(ogdf::node n1, ogdf::node n2, bool nondegenerate);

	clock_t cycles_;
};

} /* namespace mco */
#endif /* ONLINE_VERTEX_ENUMERATOR_H_ */
