//
//  ep_dual_benson.h
//  mco
//
//  Created by Fritz Bökler on 10.04.14.
//
//

#ifndef __mco__ep_dual_benson__
#define __mco__ep_dual_benson__

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/basic/abstract_solver.h>
#include <mco/basic/weight_function_adaptors.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>
#include <mco/generic/benson_dual/ove_fp_v2.h>
#include <mco/ep/basic/dijkstra.h>

namespace mco {

class LexDijkstraSolverAdaptor {
public:
    LexDijkstraSolverAdaptor(const ogdf::Graph& graph,
                             std::function<const Point *(const ogdf::edge)> weights,
                             ogdf::node source,
                             ogdf::node target)
    :   graph_(graph),
        weights_(weights),
        source_(source),
        target_(target) {}
    
    inline double operator()(const Point& weighting,
                             Point& value);

private:
    LexDijkstra lex_dijkstra_solver_;
    const ogdf::Graph& graph_;
    std::function<const Point *(const ogdf::edge)> weights_;
    const ogdf::node source_;
    const ogdf::node target_;
    
};
    
template<typename OnlineVertexEnumerator = GraphlessOVE>
class EPDualBensonSolver : public AbstractSolver {
public:
    EPDualBensonSolver(double epsilon = 1E-8)
    :   epsilon_(epsilon) {}
    
    void Solve(const ogdf::Graph& graph,
               std::function<Point const * (const ogdf::edge)> weight,
               const ogdf::node source,
               const ogdf::node target);
    
private:
    double epsilon_;
    
};
    
inline double LexDijkstraSolverAdaptor::
operator()(const Point& weighting, Point& value) {
    
    unsigned dimension = weights_(graph_.chooseEdge())->dimension();
    
    ogdf::NodeArray<Point *> distance(graph_, nullptr);
    ogdf::NodeArray<ogdf::edge> predecessor(graph_);
    
    for(auto n: graph_.nodes) {
        distance[n] = new Point(dimension + 1);
    }
    
    lex_dijkstra_solver_.singleSourceShortestPaths(graph_,
                                                   LexWeightFunctionAdaptor(graph_, weights_, weighting),
                                                   source_,
                                                   distance,
                                                   predecessor,
                                                   DijkstraModes::Undirected);
    
    Point& target_cost = *distance[target_];
    
    for(unsigned i = 0; i < dimension; ++i) {
        value[i] = target_cost[i + 1];
    }
    
    double weighted_value = target_cost[0];
    
    for(auto p : graph_.nodes) {
        delete distance[p];
    }
    
    return weighted_value;
}
 
template<typename OnlineVertexEnumerator>
inline void EPDualBensonSolver<OnlineVertexEnumerator>::
Solve(const ogdf::Graph& graph,
      std::function<Point const *(const ogdf::edge)> weights,
      ogdf::node source,
      ogdf::node target) {
    
    std::list<Point *> solutions;
    
    DualBensonScalarizer<OnlineVertexEnumerator>
    dual_benson_solver(LexDijkstraSolverAdaptor(graph, weights, source, target),
                       weights(graph.chooseEdge())->dimension(),
                       epsilon_);
    
    dual_benson_solver.Calculate_solutions(solutions);
    add_solutions(solutions.begin(), solutions.end());
    
}
    
}

#endif /* defined(__mco__ep_dual_benson__) */
