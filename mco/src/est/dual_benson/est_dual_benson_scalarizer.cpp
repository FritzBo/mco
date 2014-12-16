/*
 * est_dual_benson_scalarizer.cpp
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#include <mco/est/dual_benson/est_dual_benson_scalarizer.h>

#include <list>
#include <vector>

using std::vector;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::EdgeArray;
using ogdf::edge;
using ogdf::node;

#include <mco/basic/point.h>

namespace mco {

//template<typename OnlineVertexEnumerator>
//double ESTDualBensonScalarizer<OnlineVertexEnumerator>::
//operator()(const Point &weighting, Point &value) {
//
//    auto instance = *instance_;
//	clock_t start = clock();
//    parents_.init(instance.graph());
//	unsigned int dim = instance.dimension();
//	const EdgeArray<Point *> &weights = instance.weights();
//
//	vector<edge> sorted_edges;
//	EdgeArray<double> weighted_costs(instance.graph());
//
//	edge e;
//	forall_edges(e, *instance_().graph()) {
//		weighted_costs[e] = *weights[e] * weighting;
//		sorted_edges.push_back(e);
//	}
//
//	sort(sorted_edges.begin(), sorted_edges.end(), [this, dim, weights, weighted_costs] (edge e1, edge e2) {
//		if(weighted_costs[e1] - weighted_costs[e2] < -epsilon_)
//		    return true;
//		else if(weighted_costs[e2] - weighted_costs[e1] < -epsilon_)
//			return false;
//		else
//			for(unsigned int i = 0; i < dim; ++i) {
//				if((*weights[e1])[i] - (*weights[e2])[i] < -epsilon_)
//					return true;
//				if((*weights[e2])[i] - (*weights[e1])[i] < -epsilon_)
//					return false;
//			}
//		return false;
//	});
//
//	cycles_ += clock() - start;
//	return kruskal_solver(sorted_edges, weighted_costs, value);
//}

//template<typename OnlineVertexEnumerator>
//double ESTDualBensonScalarizer<OnlineVertexEnumerator>::
//kruskal_solver(vector<edge> &sorted_edges,
//               EdgeArray<double> costs, Point &value) {
//
//	double cost = 0;
//
//	node n;
//	forall_nodes(n, *instance_().graph()) {
//		make_set(n);
//	}
//
//	//TODO: Stop after n-1 edges have been added
//	for(edge e: sorted_edges) {
//		if(find_set(e->source()) != find_set(e->target())) {
//			set_union(e->source(), e->target());
//			cost += costs[e];
//			value += *(*instance_().weights())[e];
//		}
//	}
//
//	return cost;
//}

//template<typename OnlineVertexEnumerator>
//void ESTDualBensonScalarizer<OnlineVertexEnumerator>::make_set(ogdf::node n) {
//	parents_[n] = n;
//}

//template<typename OnlineVertexEnumerator>
//ogdf::node ESTDualBensonScalarizer<OnlineVertexEnumerator>::find_set(ogdf::node n) {
//	node current = n;
//	list<node> path_nodes;
//
//	while(parents_[current] != current) {
//		path_nodes.push_back(current);
//		current = parents_[current];
//	}
//
//	for(node n: path_nodes)
//		parents_[n] = current;
//
//	return current;
//}

//template<typename OnlineVertexEnumerator>
//void ESTDualBensonScalarizer<OnlineVertexEnumerator>::set_union(ogdf::node u, ogdf::node v) {
//	parents_[find_set(v)] = parents_[find_set(u)];
//}

} /* namespace mco */
