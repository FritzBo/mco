#pragma once
/*
 * graph_vertex_container.h
 *
 *  Created on: 14.08.2013
 *      Author: fritz
 */

#ifndef GRAPH_VERTEX_CONTAINER_H_
#define GRAPH_VERTEX_CONTAINER_H_

#include <map>
#include <set>
#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/generic/geometric/online_vertex_enumerator.h>

namespace mco {

class GraphVertexContainer: public OnlineVertexEnumerator {
public:
	GraphVertexContainer(Point &ideal_point, unsigned int dimension, double epsilon);
};

} /* namespace mco */
#endif /* GRAPH_VERTEX_CONTAINER_H_ */
