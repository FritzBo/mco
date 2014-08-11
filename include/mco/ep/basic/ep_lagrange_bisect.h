//
//  ep_lagrange_bisect.h
//  mco
//
//  Created by Fritz Bökler on 07.08.14.
//
//

#ifndef __mco__ep_lagrange_bisect__
#define __mco__ep_lagrange_bisect__

#include <functional>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class EpLagrangeBisect {
    
    using cost_function = std::function<const Point&(const ogdf::edge)>;

public:
    double find_lagrange_multi(const ogdf::Graph& graph,
                               cost_function costs,
                               unsigned dimension,
                               const ogdf::node source,
                               const ogdf::node target,
                               bool directed,
                               const Point& bounds,
                               Point& lambda);
                               

private:
    class DijkstraAdaptor {
    public:
        DijkstraAdaptor(const ogdf::Graph& graph,
                        cost_function costs,
                        unsigned dimension,
                        const ogdf::node source,
                        const ogdf::node target,
                        bool directed)
        :   graph_(graph),
            costs_(costs),
            dimension_(dimension),
            source_(source),
            target_(target),
            directed_(directed) { }
        
        double operator()(const Point& weighting,
                          Point& value);
        
    private:
        const ogdf::Graph& graph_;
        cost_function costs_;
        unsigned dimension_;
        const ogdf::node source_;
        const ogdf::node target_;
        bool directed_;
    };
};
    
}

#endif /* defined(__mco__ep_lagrange_bisect__) */
