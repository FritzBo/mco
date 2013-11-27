/*
 * online_vertex_enumerator.cpp
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#include <mco/generic/geometric/online_vertex_enumerator.h>

#include <cassert>
#include <set>
#include <list>
#include <vector>
#include <cmath>

using std::make_pair;
using std::vector;
using std::list;
using std::set;
using std::abs;

//#include <setoper.h>
//#include <cdd.h>

#include <ogdf/basic/Graph.h>

using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;
using ogdf::AdjElement;

namespace mco {

//OnlineVertexEnumerator::OnlineVertexEnumerator(std::list<Point> initial_inequalities, list<double> initial_rhs, unsigned int dimension, double epsilon = 1E-6) :
//		dimension_(dimension), epsilon_(epsilon), comp_(epsilon), point_nodes_(comp_) {
//	node_points_.init(vertex_graph_);
//
//	dd_set_global_constants();
//
//	dd_MatrixPtr M, G;
//	dd_PolyhedraPtr P;
//	dd_ErrorType err;
//
//	M = dd_CreateMatrix(dimension_, dimension_ + 1);
//
//	auto inequality_iterator = initial_inequalities.begin();
//	auto rhs_iterator = initial_rhs.begin();
//
//	unsigned int m = 0;
//	while(inequality_iterator != initial_inequalities.end()) {
//		for(unsigned int i = 0; i < dimension_; ++i)
//			dd_set_d(M->matrix[m][i + 1], (*inequality_iterator)[i]);
//		dd_set_d(M->matrix[m][0], *rhs_iterator);
//
//		m++;
//		inequality_iterator++;
//		rhs_iterator++;
//		if(m >= dimension_)
//			break;
//	}
//
//	M->representation = dd_Inequality;
//	P = dd_DDMatrix2Poly(M, &err);
//	G = dd_CopyGenerators(P);
//
//	node finite_node, n;
//	Point *finite_point;
//	list<node> infinity_points;
//
//	cout << G->rowsize << endl;
//
//	for(unsigned int i = 0; i < G->rowsize; ++i) {
//		Point *p = new Point(dimension_ + 1);
//		for(unsigned int j = 1; j < dimension_ + 1; ++j)
//			(*p)[j - 1] = dd_get_d(G->matrix[i][j]);
//		(*p)[dimension_] = dd_get_d(G->matrix[i][0]);
//
//		cout << *p << endl;
//
//		n = vertex_graph_.newNode();
//
//		if((*p)[dimension_] == 1) {
//			finite_node = n;
//			finite_point = p;
//		} else
//			infinity_points.push_back(n);
//
//		node_points_[n] = p;
//		point_nodes_.insert(make_pair(p, n));
//	}
//
//	for(node n: infinity_points) {
//		vertex_graph_.newEdge(n, finite_node);
//		vertex_graph_.newEdge(finite_node, n);
//	}
//
//	unprocessed_projective_points_.push_back(finite_point);
//
//	dd_FreeMatrix(M);
//	dd_FreeMatrix(G);
//	dd_FreePolyhedra(P);
//
//	dd_free_global_constants();
//
//	while(inequality_iterator != initial_inequalities.end()) {
//		node n;
//		Point *p;
//		bool broken = false;
//		forall_nodes(n, vertex_graph_) {
//			p = node_points_[n];
//			if((*inequality_iterator) * (*p) - (*rhs_iterator) * (*p)[dimension_] < -epsilon) {
//				broken = true;
//				break;
//			}
//		}
//
//		if(broken)
//			add_hyperplane(*p, *inequality_iterator, *rhs_iterator);
//
//		inequality_iterator++;
//		rhs_iterator++;
//	}
//}

OnlineVertexEnumerator::~OnlineVertexEnumerator() {
	node n;

	forall_nodes(n, vertex_graph_)
		delete node_points_[n];

	for(auto inequality : list_of_inequalities_)
		delete inequality;

}

Point* OnlineVertexEnumerator::to_projective(Point &non_projective_point) {
	double * projective_values = new double[dimension_ + 1];

	for(unsigned int i = 0; i < dimension_; ++i)
		projective_values[i] = non_projective_point[i];

	projective_values[dimension_] = 1;

	Point *projective_point = new Point(projective_values, dimension_ + 1);
	delete[] projective_values;

	return projective_point;
}

Point OnlineVertexEnumerator::normalize_projective(Point projective_point) {
	double *normalized_projective_values = new double[dimension_ + 1];

	for(unsigned int i = 0; i < dimension_; ++i)
		normalized_projective_values[i] = projective_point[i] / projective_point[dimension_];

	normalized_projective_values[dimension_] = 1;

	Point normalized_projective_point(normalized_projective_values, dimension_ + 1);
	delete[] normalized_projective_values;

	return normalized_projective_point;
}

node OnlineVertexEnumerator::get_node(Point &non_projective_point) {
	Point *projective_point = to_projective(non_projective_point);
	node n = point_nodes_[projective_point];
//	cout << n << endl;

//	for(auto entity: point_nodes_)
//		cout << *entity.first << " - " << entity.second << endl;

	assert(node_points_[n]->is_equal(*projective_point, epsilon_));

	delete projective_point;


	return n;
}

bool OnlineVertexEnumerator::has_next() {
	while(!unprocessed_projective_points_.empty() && point_nodes_.count(unprocessed_projective_points_.front()) == 0)
			unprocessed_projective_points_.pop_front();

	return !unprocessed_projective_points_.empty();
}

Point * OnlineVertexEnumerator::next_vertex() {
	Point *point;

//	cout << *unprocessed_projective_points_.front() << " noch im Graphen: " << point_nodes_.count(unprocessed_projective_points_.front()) << endl;
	while(!unprocessed_projective_points_.empty() && point_nodes_.count(unprocessed_projective_points_.front()) == 0)
		unprocessed_projective_points_.pop_front();

//	cout << *node_points_[point_nodes_[unprocessed_projective_points_.front()]] << endl;

	if(unprocessed_projective_points_.empty())
		return nullptr;

	point = unprocessed_projective_points_.front();
//	cout << *point << ", dim: " << point->dimension() << endl;
	Point * new_point = new Point(dimension_);

//	for(auto entity: point_nodes_)
//		cout << *entity.first << " - " << entity.second << ", dim: " << entity.first->dimension() << endl;

	assert(point->dimension() == dimension_ + 1);

	if(abs((*point)[dimension_] - 1) > epsilon_) {
		cout << "Point " << *point << " is not normalized: " << (*point)[dimension_] - 1 << endl;
		assert(false);
	}

	for(unsigned int i = 0; i < dimension_; ++i)
		(*new_point)[i] = (*point)[i];

	unprocessed_projective_points_.pop_front();

	return new_point;
}

void OnlineVertexEnumerator::add_hyperplane(Point &vertex, Point &normal, double rhs) {
//	cout << vertex << endl;
//	cout << normal << endl;
//	cout << rhs << endl;

	clock_t start = clock();

	Point *projective_normal = to_projective(normal);
	(*projective_normal)[dimension_] = -rhs;

	list_of_inequalities_.push_back(projective_normal);

	// Check if redundancy is introduced
//	dd_set_global_constants();
//	dd_MatrixPtr A;
//	dd_ErrorType err;
//	dd_rowset rows;
//	A = dd_CreateMatrix(list_of_inequalities_.size(), dimension_ + 1);
//	int j = 0;
//	for(auto row: list_of_inequalities_) {
//
//		for(unsigned int i = 0; i < dimension_; ++i)
//			dd_set_d(A->matrix[j][i + 1], (*row)[i]);
//
//		dd_set_d(A->matrix[j][0], (*row)[dimension_]);
//
//		j+= 1;
//	}
//	A->representation = dd_Inequality;
//	rows = dd_RedundantRows(A, &err);
//
//	if(set_card(rows) > 0) {
//		cout << "CDD giving redundant inequalities: " << set_card(rows) << endl;
//		dd_WriteMatrix(stdout, A);
//		set_write(rows);
//		assert(false);
//	}
//	set_free(rows);
//	dd_FreeMatrix(A);
//	dd_free_global_constants();
	// End Redundancy check


	list<node> active_nodes;
	active_nodes.push_back(get_node(vertex));

	list<node> new_face_nodes;

	NodeArray<bool> already_active(vertex_graph_, false);

	bool nondegenerate = true;

	while(!active_nodes.empty()) {

		node active_node = active_nodes.front();
		active_nodes.pop_front();
//		cout << "active node: " << active_node << endl;
		assert(active_node != nullptr);

		Point *projective_vertex = node_points_[active_node];

//		cout << "vertex: " << *projective_vertex << endl;

		NodeArray<bool> neighor_checked(vertex_graph_, false);

		AdjElement *adj;
		forall_adj(adj, active_node) {

			node neighbor = adj->theEdge()->target();

//			cout << "neighbor: " << neighbor << endl;

			if(neighbor == active_node)
				continue;

			if(neighor_checked[neighbor])
				continue;

			neighor_checked[neighbor] = true;

			Point *projective_check_point = node_points_[neighbor];

//			cout << "neighbor to check: " << *projective_check_point << endl;

//			for(auto n : new_face_nodes)
//				cout << "1 face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

			node new_node = nullptr;

			// Check if neighbor is on the other side of the hyperplane
			if((*projective_normal) * (*projective_check_point) < -epsilon_) {
				if(!already_active[neighbor]) {
					active_nodes.push_back(neighbor);
					already_active[neighbor] = true;
				}
//				cout << "neighbor is also outside" << endl;
				continue;

			} else if((*projective_normal) * (*projective_check_point) < epsilon_) {
				nondegenerate = false;
				new_node = neighbor;
//				cout << "neighbor " << new_node << " resides on hyperplane" << endl;
				bool already_on_face = false;
				for(auto n : new_face_nodes) {
					if(new_node->index() == n->index())
						already_on_face = true;
				}
				if(already_on_face)
					continue;

			} else {

//				cout << "On hyperplane check: " << (*projective_normal) * (*projective_check_point) << endl;

				double alpha = - (*projective_normal) * (*projective_vertex);
				Point diff_direction = (*projective_check_point) - (*projective_vertex);

				if(abs((*projective_normal) * diff_direction) < epsilon_) {
//					cout << "divisor: " << (*projective_normal) * diff_direction << endl;
					continue;
				} else
					alpha /=  (*projective_normal) * diff_direction;

//				for(auto n : new_face_nodes)
//					cout << "2 face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

//				cout << "diff direction: " << diff_direction << endl;
//				cout << "alpha: " << alpha << endl;

				assert( alpha > epsilon_ && alpha < 1 + epsilon_);

//				for(auto n : new_face_nodes)
//					cout << "3 face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

				diff_direction *= alpha;
				Point *projective_cut_point = new Point(normalize_projective((*projective_vertex) + diff_direction));

//				cout << "new cut vertex: " << *projective_cut_point << endl;

				unprocessed_projective_points_.push_back(projective_cut_point);

//				for(auto n : new_face_nodes)
//					cout << "4 face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

//				cout << adj->theEdge() << endl;
				edge new_edge = vertex_graph_.split(adj->theEdge());
//				cout << adj->theEdge() << endl;
				new_node = new_edge->source();
//				cout << (void *) new_node << endl;

//				cout << "New node: " << new_node << endl;

				vertex_graph_.newEdge(neighbor, new_node);

				node_points_[new_node] = projective_cut_point;

//				for(auto n : new_face_nodes)
//					cout << "5 face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

				point_nodes_.insert(make_pair(projective_cut_point, new_node));
			}

//			for(auto n : new_face_nodes)
//				cout << "face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

			new_face_nodes.push_back(new_node);

//			for(auto n : new_face_nodes)
//				cout << "face node: " << n << " with point: " << *node_points_[n] << " at adress " << (void *) n << endl;

		}	// forall_adj

//		cout << "Deleting node: " << active_node << endl;
		vertex_graph_.delNode(active_node);
		point_nodes_.erase(projective_vertex);

	}	// active node

//	cout << "Checking face nodes: ";
//	for(auto n : new_face_nodes)
//		cout << n << ", ";
//	cout << endl;

	for(auto n : new_face_nodes) {
		for(auto node : new_face_nodes) {

			if(n == node)
				continue;

			if(inside_face(*node_points_[n], *node_points_[node], nondegenerate))
				continue;

			vertex_graph_.newEdge(node, n);
			vertex_graph_.newEdge(n, node);
		}

		if(dimension_ == 2)
			assert((n->indeg() >= 2 && n->outdeg() >= 2) || abs((*node_points_[n])[dimension_]) < epsilon_);
		else
			if((n->indeg() < 3 || n->outdeg() < 3) && abs((*node_points_[n])[dimension_]) > epsilon_) {
				cout << "node " << n << " with point " << *node_points_[n] << " has indegree " << n->indeg() << " and outdegree " << n->outdeg() << endl;
				assert(false);
			}
	}

	cycles_ += clock() - start;
}

unsigned int OnlineVertexEnumerator::number_of_hyperplanes() {
	return list_of_inequalities_.size();
}

bool OnlineVertexEnumerator::inside_face(Point& p1, Point& p2, bool nondegenerate) {
//	cout << "p1: " << p1 << " (" << point_nodes_[&p1] << "), p2: " << p2 << " (" << point_nodes_[&p2] << ")" << endl;
	if(dimension_ == 2)
		return false;

	else if(dimension_ > 3 || !nondegenerate) {

		assert(!p1.is_equal(p2, epsilon_));

		if(abs(p1[dimension_]) < epsilon_ && abs(p2[dimension_]) < epsilon_)
			return false;

		set<Point *, LexicographicPointComparator> common_vertices(comp_);
		unsigned int i = 0;
		for(auto inequality : list_of_inequalities_) {

			set<Point *, LexicographicPointComparator> new_vertices(comp_);

//			cout << "inequality: " << *inequality << endl;

			if(		abs((*inequality) * p1) < epsilon_ &&
					abs((*inequality) * p2) < epsilon_) {

				node n;
				forall_nodes(n, vertex_graph_) {
//					cout << "checking node " << n << " with point " << *node_points_[n] << ": " << (*inequality) * (*node_points_[n]) << endl;
					if(abs((*inequality) * (*node_points_[n])) < epsilon_)
						new_vertices.insert(node_points_[n]);
				}

				if(common_vertices.empty())
					common_vertices.insert(new_vertices.begin(), new_vertices.end());
				else {
					list<Point *> temp_points;
					set_intersection(common_vertices.begin(), common_vertices.end(), new_vertices.begin(), new_vertices.end(), back_inserter(temp_points), comp_);
//					cout << "intersection size: "<< temp_points.size() << endl;
					common_vertices.clear();
					common_vertices.insert(temp_points.begin(), temp_points.end());
				}

//				cout << "current number of common vertices: " << common_vertices.size() << endl;

			} else
				assert(i < list_of_inequalities_.size() - 1);

			i += 1;
		}

//		cout << "Number of common vertices: " << common_vertices.size() << endl;

		if(common_vertices.size() == 2)
			return false;
		else
			return true;

	} else if(dimension_ == 3 || nondegenerate) {

		unsigned int tight_inequalities = 0;
		for(auto inequality : list_of_inequalities_) {
//			cout << abs((*inequality) * p1) << " && " << abs((*inequality) * p2) << endl;
			if(abs((*inequality) * p1) < epsilon_ && abs((*inequality) * p2) < epsilon_)
				tight_inequalities++;
		}

//		cout << "tight inequalities: " << tight_inequalities << endl;

		if(tight_inequalities == dimension_ - 1)
			return false;
		else if(tight_inequalities < dimension_ - 1)
			return true;
		else
			assert(false);
	}

	return false;
}

} /* namespace mco */
