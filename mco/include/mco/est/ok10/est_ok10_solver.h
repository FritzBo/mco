//
//  est_ok10_solver.h
//  mco
//
//  Created by Fritz Bökler on 20.1.2015
//
//

#ifndef mco_est_ok10_solver_h
#define mco_est_ok10_solver_h

#include <functional>

#include <mco/basic/abstract_solver.h>
#include <mco/basic/weight_function_adaptors.h>
#include <mco/est/basic/kruskal_st_solver.h>
#include <mco/generic/ok10/ok10_scalarizer_projective.h>

namespace mco {

class LexKruskalSolverAdaptor {
public:
    inline LexKruskalSolverAdaptor(ogdf::Graph& graph,
                                   ogdf::EdgeArray<Point *>& costs,
                                   unsigned objectives);

    inline double operator()(const Point& weighting, Point& value);

private:
    LexKruskalSTSolver lex_est_solver_;
    ogdf::Graph& graph_;
    std::function<const Point*(const ogdf::edge)> cost_function_;
    unsigned objectives_;
};


class EstOk10Solver
: public AbstractSolver<std::list<ogdf::edge>> {

public:
    EstOk10Solver(double epsilon = 1E-8)
    :   epsilon_(epsilon) { }

    void Solve(Graph& graph,
               EdgeArray<Point *>& costs,
               unsigned objectives);

    unsigned number_facets() {
        return facets_;
    }

private:
    double epsilon_;
    unsigned facets_;
};

LexKruskalSolverAdaptor::
LexKruskalSolverAdaptor(ogdf::Graph& graph,
                        ogdf::EdgeArray<Point *>& costs,
                        unsigned objectives)
:   graph_(graph),
    objectives_(objectives) {

    cost_function_ = [costs] (ogdf::edge e) {
        return costs(e);
    };
}

inline double LexKruskalSolverAdaptor::
operator()(const Point& weighting,
           Point& value) {

    lex_est_solver_.Solve(graph_,
            LexWeightFunctionAdaptor(graph_,
                                     cost_function_,
                                     weighting),
                          objectives_ + 1);

    Point result = lex_est_solver_.min_cost();
    
    for(unsigned i = 0; i < objectives_; ++i) {
        value[i] = result[i + 1];
    }
    
    return result[0];
}
    
inline void EstOk10Solver::
Solve(Graph& graph,
      EdgeArray<Point *>& weights,
      unsigned objectives) {

    std::list<Point *> frontier;

    Point upper_bounds(-std::numeric_limits<double>::infinity(), objectives);

    for(auto e : graph.edges) {
        for(unsigned i = 0; i < objectives; ++i) {
            upper_bounds[i] = max(upper_bounds[i],
                                  weights(e)->operator[](i));
        }
    }

    for(unsigned i = 0; i < objectives; ++i) {
        upper_bounds[i] *= graph.numberOfNodes() / 2;
    }

    LexKruskalSolverAdaptor solver(graph, weights, objectives);

    Ok10ScalarizerProjective scalarizer(solver,
                                        upper_bounds,
                                        objectives);

    scalarizer.Calculate_solutions(frontier);

    std::list<std::pair<std::list<ogdf::edge>, Point>> solutions;

    for(auto point : frontier) {
        solutions.push_back(make_pair(std::list<ogdf::edge>(), *point));
    }

    add_solutions(solutions.begin(), solutions.end());
}
    
    
} /* namespace mco */

#endif
