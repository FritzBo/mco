//
//  ep_lagrange_bisect.cpp
//  mco
//
//  Created by Fritz Bökler on 07.08.14.
//
//

#include <mco/ep/basic/ep_lagrange_bisect.h>

using ogdf::edge;
using ogdf::node;

#include <mco/generic/lagrange_bisect.h>

#include <mco/ep/basic/dijkstra.h>

namespace mco {

double EpLagrangeBisect::find_lagrange_multi(const ogdf::Graph& graph,
                                             cost_function costs,
                                             unsigned dimension,
                                             const ogdf::node source,
                                             const ogdf::node target,
                                             bool directed,
                                             const Point& bounds,
                                             Point& lambda) {
    
    LagrangeBisect bisect;
    
    return bisect.find_multiplier(DijkstraAdaptor(graph,
                                                  costs,
                                                  dimension,
                                                  source,
                                                  target,
                                                  directed),
                                  dimension,
                                  bounds,
                                  lambda);
}

double EpLagrangeBisect::DijkstraAdaptor::
operator()(const Point& weighting,
           Point& value) {
    
    ogdf::NodeArray<double> distance(graph_);
    ogdf::NodeArray<ogdf::edge> predecessor(graph_);
    
    Dijkstra<double, PairComparator<double, std::less<double>>> sssp_solver;
    
    auto weight = [&weighting, this] (edge e) {
        return this->costs_(e) * weighting;
    };
    
    sssp_solver.singleSourceShortestPaths(graph_,
                                          weight,
                                          source_,
                                          predecessor,
                                          distance,
                                          directed_ ? DijkstraModes::Forward : DijkstraModes::Undirected);
    
    for(unsigned i = 0; i < dimension_; ++i) {
        value[i] = 0;
    }
    
    node current_node = target_;
    edge predecessor_edge;
    while(current_node != source_) {
        predecessor_edge = predecessor[current_node];
        value += costs_(predecessor_edge);
        current_node = current_node == predecessor_edge->target() ? predecessor_edge->source() : predecessor_edge->target();
    }
                                          
    return distance[target_];
    
}
    
}