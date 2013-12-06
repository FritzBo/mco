/*
 * dual_benson_vertex_container.cpp
 *
 *  Created on: 01.10.2013
 *      Author: fritz
 */

#include <mco/generic/benson_weightspace/dual_benson_vertex_container.h>

#include <list>

using std::list;
using std::make_pair;

#include <ogdf/basic/Graph.h>

using ogdf::node;

#include <mco/point.h>

namespace mco {

DualBensonVertexContainer::DualBensonVertexContainer(Point &initial_value, unsigned int dimension, double epsilon) :
	OnlineVertexEnumerator(dimension, epsilon) {

	node n, v;
	Point *p;
	for(unsigned int i = 0; i < dimension_ - 1; ++i) {
		n = vertex_graph_.newNode();
		node_inequality_indices_[n] = new list<int>();
		p = new Point(dimension_ + 1);
		for(unsigned int j = 0; j < dimension_ - 1; ++j) {
			(*p)[j] = i == j ? 1 : 0;

			if(i != j)
				node_inequality_indices_[n]->push_back(j);
		}
		(*p)[dimension_ - 1] = initial_value[i];
		(*p)[dimension_] = 1;

		node_inequality_indices_[n]->push_back(dimension_ - 1);
		node_inequality_indices_[n]->push_back(dimension_);
		birth_index_[n] = dimension_;

		point_nodes_.insert(make_pair(p, n));
		node_points_[n] = p;
		unprocessed_projective_points_.push(p);

		p = new Point(dimension_ + 1);
		for(unsigned int j = 0; j < dimension_; ++j)
			(*p)[j] = i == j ? 1 : 0;
		(*p)[dimension_] = 0;

		list_of_inequalities_.push_back(p);

		forall_nodes(v, vertex_graph_) {
			if(v != n) {
				vertex_graph_.newEdge(v, n);
				vertex_graph_.newEdge(n, v);
			}
		}
	}

	p = new Point(dimension_ + 1);
	for(unsigned int j = 0; j < dimension_ - 1; ++j)
		(*p)[j] = -1;
	(*p)[dimension_ - 1] = 0;
	(*p)[dimension_] = 1;

	list_of_inequalities_.push_back(p);

	n = vertex_graph_.newNode();
	node_inequality_indices_[n] = new list<int>();
	p = new Point(dimension_ + 1);
	for(unsigned int j = 0; j < dimension_; ++j)
		(*p)[j] = 0;

	(*p)[dimension_ - 1] = initial_value[dimension_ - 1];
	(*p)[dimension_] = 1;

	for(unsigned int i = 0; i < dimension_ - 1; ++i)
		node_inequality_indices_[n]->push_back(i);

	node_inequality_indices_[n]->push_back(dimension_);
	birth_index_[n] = dimension_;

	point_nodes_.insert(make_pair(p, n));
	node_points_[n] = p;
	unprocessed_projective_points_.push(p);

	forall_nodes(v, vertex_graph_) {
		if(v != n) {
			vertex_graph_.newEdge(v, n);
			vertex_graph_.newEdge(n, v);
		}
	}

	p = new Point(dimension_ + 1);
	for(unsigned int j = 0; j < dimension_; ++j)
		(*p)[j] = initial_value[j] - initial_value[dimension_ - 1];
	(*p)[dimension_ - 1] = -1;
	(*p)[dimension_] = initial_value[dimension_ - 1];

	list_of_inequalities_.push_back(p);

	n = vertex_graph_.newNode();
	node_inequality_indices_[n] = new list<int>();
	p = new Point(dimension_ + 1);
	for(unsigned int j = 0; j < dimension_ - 1; ++j)
		(*p)[j] = 0;
	(*p)[dimension_ - 1] = -1;
	(*p)[dimension_] = 0;

	for(unsigned int i = 0; i < dimension_; ++i)
		node_inequality_indices_[n]->push_back(i);

	birth_index_[n] = dimension_ - 1;

	point_nodes_.insert(make_pair(p, n));
	node_points_[n] = p;

	forall_nodes(v, vertex_graph_) {
		if(v != n) {
			vertex_graph_.newEdge(v, n);
			vertex_graph_.newEdge(n, v);
		}
	}

//	cout << "inequalities:" << endl;
//	for(auto ineq : list_of_inequalities_)
//		cout << *ineq << endl;
//
//	cout << "vertices:" << endl;
//	forall_nodes(v, vertex_graph_) {
//		cout << *node_points_[v] << endl;
//		cout << "Node inequalities:" << endl;
//			for(auto index: *node_inequality_indices_[v])
//				cout << index << ", ";
//			cout << endl;
//	}
//
//	cout << "graph has " << vertex_graph_.numberOfNodes() << " nodes and " << vertex_graph_.numberOfEdges() << " edges" << endl;
}

DualBensonVertexContainer::~DualBensonVertexContainer() {
}

} /* namespace mco */
