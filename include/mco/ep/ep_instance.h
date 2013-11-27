#pragma once
/*
 * EpInstance.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef EPINSTANCE_H_
#define EPINSTANCE_H_

#include <exception>
#include <string>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/EdgeArray.h>

#include "../vector_weighted_graph.h"

namespace mco {

class Point;

class EpInstance : public VectorWeightedGraph {

public:
	EpInstance(std::shared_ptr<ogdf::Graph> graph, std::shared_ptr<ogdf::EdgeArray<Point *>> weights, unsigned int dimension, ogdf::node const source, ogdf::node const target);

	ogdf::NodeElement * source() const { return source_; }
	ogdf::NodeElement * target() const { return target_; }

private:
	ogdf::NodeElement * const source_;
	ogdf::NodeElement * const target_;
};

}

#endif /* EPINSTANCE_H_ */
