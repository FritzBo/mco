#pragma once
/*
 * vector_weighted_graph.h
 *
 *  Created on: 25.03.2013
 *      Author: fritz
 */

#ifndef VECTOR_WEIGHTED_GRAPH_H_
#define VECTOR_WEIGHTED_GRAPH_H_

#include <memory>

#include <ogdf/basic/Graph.h>

#include "point.h"

namespace mco {

class VectorWeightedGraph {
	std::shared_ptr<ogdf::Graph> graph_;
	std::shared_ptr<ogdf::EdgeArray<Point *>> weights_;
	unsigned int dimension_;

public:
	VectorWeightedGraph(std::shared_ptr<ogdf::Graph> graph, std::shared_ptr<ogdf::EdgeArray<Point *>> weights, unsigned int dimension);

	VectorWeightedGraph(const VectorWeightedGraph & graph) = delete;
	void operator=(const VectorWeightedGraph & graph) = delete;

	unsigned int dimension() const { return dimension_; }
	std::shared_ptr<ogdf::EdgeArray<Point *>>  weights() { return weights_; }
	std::shared_ptr<ogdf::Graph> graph() { return graph_; }
};

} /* namespace mco */

#endif /* VECTOR_WEIGHTED_GRAPH_H_ */
