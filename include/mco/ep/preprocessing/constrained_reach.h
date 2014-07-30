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
                    unsigned dimension,
                    std::function<double(ogdf::node, unsigned)> heuristic,
                    std::list<Point> linear_bounds);
private:
};
    
}

#endif /* defined(__mco__constrained_reach__) */
