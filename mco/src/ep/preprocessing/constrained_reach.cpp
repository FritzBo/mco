//
//  constrained_reach.cpp
//  mco
//
//  Created by Fritz Bökler on 30.07.14.
//
//

#include <mco/ep/preprocessing/constrained_reach.h>

using std::function;
using std::list;
using std::vector;

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;

#include <mco/ep/basic/dijkstra.h>

using mco::Dijkstra;
using mco::DijkstraModes;

namespace mco {

void ConstrainedReachabilityPreprocessing::preprocess(Graph &graph,
                                                      function<Point*(ogdf::edge)> cost_function,
                                                      unsigned dimension,
                                                      const ogdf::node source,
                                                      const ogdf::node target,
                                                      list<Point> linear_bound,
                                                      bool directed) {
    
    list<ogdf::node> candidates;
    
    Dijkstra<double, PairComparator<double, std::less<double>>> sssp_solver;
    
    NodeArray<ogdf::edge> predecessor(graph);
    vector<NodeArray<double>> source_distances(dimension, graph);
    vector<NodeArray<double>> target_distances(dimension, graph);
    
    for(unsigned i = 0; i < dimension; ++i) {
        
        
        auto weight = [cost_function, i] (ogdf::edge e) {return (cost_function(e))->operator[](i); };
    
        sssp_solver.singleSourceShortestPaths(graph,
                                              weight,
                                              source,
                                              predecessor,
                                              source_distances[i],
                                              directed ? DijkstraModes<>::Forward : DijkstraModes<>::Undirected);
        
        sssp_solver.singleSourceShortestPaths(graph,
                                              weight,
                                              target,
                                              predecessor,
                                              target_distances[i],
                                              directed ? DijkstraModes<>::Backward : DijkstraModes<>::Undirected);

    }

    
    for(auto n : graph.nodes) {
        for(auto bound : linear_bound) {
            double inner_product = 0;
            for(unsigned i = 0; i < dimension; ++i) {
                inner_product += bound[i] * (source_distances[i][n] + target_distances[i][n]);
            }
            
            if(inner_product > -bound[dimension]) {
                candidates.push_back(n);
                break;
            }
            
        }
    }
    
    unsigned nodes_deleted = 0;
    for(auto n : candidates) {
        graph.delNode(n);
        nodes_deleted++;
    }

#ifndef NDEBUG
    std::cout << "Numer of nodes deleted: " << nodes_deleted << std::endl;
#endif
}

}