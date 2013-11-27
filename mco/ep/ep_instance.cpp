/*
 * EpInstance.cpp
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#include <mco/ep/ep_instance.h>

#include <memory>

using std::shared_ptr;

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EdgeArray.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::EdgeArray;

using mco::EpInstance;

namespace {

bool nodeExists(shared_ptr<const Graph> graph, node no) {
	node n;
	bool exists = false;
	forall_nodes(n, *graph) {
		if(n == no) {
			exists = true;
			break;
		}
	}
	return exists;
}

}

EpInstance::EpInstance(shared_ptr<Graph> graph, shared_ptr<EdgeArray<Point *>> weights, unsigned int dimension, node const source, node const target) : VectorWeightedGraph(graph, weights, dimension), source_(source), target_(target) {
	// check source and target nodes
	if(source == nullptr)
		throw "Source node must not be NULL.";

	if(!nodeExists(graph, source))
		throw "Node does not exists!";

	if(target == nullptr)
		throw "Target node must not be NULL.";

	if(!nodeExists(graph, target))
		throw "Node does not exists!";

	// check weights
	if(!weights->valid())
		throw "Weights must be valid!";

	if(weights->graphOf() != &*graph)
		throw "Weights must be of the given graph.";

}
