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

using ogdf::Graph;
using ogdf::node;

namespace mco {

void ConstrainedReachabilityPreprocessing::preprocess(Graph &graph,
                                                      unsigned dimension,
                                                      function<double(node, unsigned)> heuristic,
                                                      list<Point> linear_bound) {
    
    list<node> candidates;
    
    for(auto n : graph.nodes) {
        for(auto bound : linear_bound) {
            double inner_product = 0;
            for(unsigned i = 0; i < dimension; ++i) {
                inner_product += bound[i] * heuristic(n, i);
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
    
    std::cout << "Numer of nodes deleted: " << nodes_deleted << std::endl;
}

}