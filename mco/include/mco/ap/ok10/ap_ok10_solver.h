//
//  ap_ok10_solver.h
//  mco
//
//  Created by Fritz Bökler on 16.11.14.
//
//

#ifndef mco_ap_ok10_solver_h
#define mco_ap_ok10_solver_h

#include <functional>

#include <mco/basic/abstract_solver.h>
#include <mco/basic/weight_function_adaptors.h>
#include <mco/ap/basic/abstract_ap_solver.h>
#include <mco/ap/basic/ap_instance.h>
#include <mco/ap/basic/lex_hungarian.h>
#include <mco/generic/ok10/ok10_scalarizer.h>

namespace mco {

class LexHungarianSolverAdaptor {
public:
    inline LexHungarianSolverAdaptor(AssignmentInstance& ap_instance);

    inline double operator()(const Point& weighting, Point& value);

private:
    LexHungarianMethod lex_ap_solver_;
    AssignmentInstance& ap_instance_;
};


class ApOk10Solver
: public AbstractSolver<std::list<ogdf::edge>> {

public:
    ApOk10Solver(double epsilon = 1E-8)
    :   epsilon_(epsilon) { }

    void Solve(AssignmentInstance & instance);

    unsigned number_facets() {
        return facets_;
    }

private:
    double epsilon_;
    unsigned facets_;
};

inline LexHungarianSolverAdaptor::
LexHungarianSolverAdaptor(AssignmentInstance& ap_instance)
:   ap_instance_(ap_instance) {
}


inline double LexHungarianSolverAdaptor::
operator()(const Point& weighting, Point& value) {

    Point result = lex_ap_solver_.Solve(ap_instance_.graph(),
                                        LexWeightFunctionAdaptor(ap_instance_.graph(),
                                                                 ap_instance_.weights(),
                                                                 weighting),
                                        ap_instance_.dimension() + 1,
                                        ap_instance_.agents());
    
    for(unsigned i = 0; i < ap_instance_.dimension(); ++i) {
        value[i] = result[i + 1];
    }
    
    return result[0];
}
    
inline void ApOk10Solver::
Solve(AssignmentInstance & instance) {

    unsigned dimension = instance.dimension();

    std::list<Point *> frontier;

    Point upper_bounds(-std::numeric_limits<double>::infinity(), dimension);

    for(auto e : instance.graph().edges) {
        for(unsigned i = 0; i < dimension; ++i) {
            upper_bounds[i] = max(upper_bounds[i],
                                  instance.weights()(e)->operator[](i));
        }
    }

    for(unsigned i = 0; i < dimension; ++i) {
        upper_bounds[i] *= instance.graph().numberOfNodes() / 2;
    }

    LexHungarianSolverAdaptor solver(instance);

    Ok10Scalarizer scalarizer(solver,
                              upper_bounds,
                              dimension);

    scalarizer.Calculate_solutions(frontier);

    std::list<std::pair<std::list<ogdf::edge>, Point>> solutions;

    for(auto point : frontier) {
        solutions.push_back(make_pair(std::list<ogdf::edge>(), *point));
    }

    add_solutions(solutions.begin(), solutions.end());
}
    
    
} /* namespace mco */

#endif
