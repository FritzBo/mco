//
//  constrained_reach.h
//  mco
//
//  Created by Fritz Bökler on 30.07.14.
//
//

#ifndef __mco__constrained_reach__
#define __mco__constrained_reach__

#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class ConstrainedReachabilityPreprocessing {
public:
    
    void preprocess(ogdf::Graph& graph,
                    std::function<Point*(ogdf::edge)> cost_function,
                    unsigned dimension,
                    const ogdf::node source,
                    const ogdf::node target,
                    std::list<Point> linear_bounds,
                    bool directed = false);
    
    void preprocess(ogdf::Graph& graph,
                    std::function<Point&(ogdf::edge)> cost_function,
                    unsigned dimension,
                    const ogdf::node source,
                    const ogdf::node target,
                    std::list<Point> linear_bounds,
                    bool directed = false) {
        
        auto cost_function_adaptor = [cost_function] (ogdf::edge e) {
            return &cost_function(e);
        };
        
        preprocess(graph,
                   cost_function_adaptor,
                   dimension,
                   source,
                   target,
                   linear_bounds);
    }
private:
};
    
}

#endif /* defined(__mco__constrained_reach__) */
