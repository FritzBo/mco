#pragma once
/*
 * assignment_instance.h
 *
 *  Created on: 19.08.2013
 *      Author: fritz
 */

#ifndef ASSIGNMENT_INSTANCE_H_
#define ASSIGNMENT_INSTANCE_H_

#include <set>
#include <memory>

#include <ogdf/basic/Graph.h>

#include <mco/vector_weighted_graph.h>

namespace mco {

class AssignmentInstance : public VectorWeightedGraph {

public:
	AssignmentInstance(std::shared_ptr<ogdf::Graph> graph, std::shared_ptr<ogdf::EdgeArray<Point *>> weights, std::shared_ptr<std::set<ogdf::node>> agents, unsigned int dimension) :
		VectorWeightedGraph(graph, weights, dimension), agents_(agents) { }

	std::shared_ptr<std::set<ogdf::node>> agents() {
		return agents_;
	}

private:
	std::shared_ptr<std::set<ogdf::node>> agents_;
};

} /* namespace mco */
#endif /* ASSIGNMENT_INSTANCE_H_ */
