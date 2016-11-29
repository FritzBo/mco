//
//  instance_scalarizer.cpp
//  mco
//
//  Created by Fritz Bökler on 29.08.14.
//
//

#include <mco/ep/basic/instance_scalarizer.h>

using ogdf::Graph;
using ogdf::EdgeArray;

namespace mco {

void InstanceScalarizer::
scaleround_instance(const Graph& graph,
                    const EdgeArray<Point>& cost_array,
                    unsigned dimension,
                    const Point& factor,
                    EdgeArray<Point>& new_cost_array,
                    bool round) {

    for(ogdf::edge edge : graph.edges) {
        for(unsigned i = 0; i < dimension; ++i) {
            if(round) {
                new_cost_array[edge][i] = std::round(cost_array[edge][i] * factor[i]);
            } else {
                new_cost_array[edge][i] = cost_array[edge][i] * factor[i];
            }
        }
    }
}

void InstanceScalarizer::
scaleround_instance(const ForwardStar& graph,
                    const FSEdgeArray<Point>& cost_array,
                    unsigned dimension,
                    const Point& factor,
                    FSEdgeArray<Point>& new_cost_array,
                    bool round) {

    for(mco::edge edge : graph.edges) {
        for(unsigned i = 0; i < dimension; ++i) {
            if(round) {
                new_cost_array[edge][i] = std::round(cost_array[edge][i] * factor[i]);
            } else {
                new_cost_array[edge][i] = cost_array[edge][i] * factor[i];
            }
        }
    }
}


}