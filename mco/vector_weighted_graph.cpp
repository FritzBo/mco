/*
 * vector_weighted_graph.cpp
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#include <mco/vector_weighted_graph.h>

#include <memory>

using std::shared_ptr;

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EdgeArray.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::EdgeArray;


namespace mco {

VectorWeightedGraph::VectorWeightedGraph(shared_ptr<Graph> graph, shared_ptr<EdgeArray<Point *>> weights, unsigned int dimension) :
		graph_(graph), weights_(weights), dimension_(dimension) {
	// check weights
	if(!weights->valid())
		throw "Weights must be valid!";

	if(weights->graphOf() != &*graph)
		throw "Weights must be of the given graph.";
}

} /* namespace mco */
