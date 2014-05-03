#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef MARTINS_B_H_
#define MARTINS_B_H_

#include <mco/basic/abstract_solver.h>

namespace mco {

class EpSolverMartins : public AbstractSolver {

public:
	explicit EpSolverMartins(double epsilon = 0)
    : epsilon_(epsilon) { }
    
	virtual void Solve(ogdf::Graph& graph,
                       std::function<const Point*(ogdf::edge)> weights,
                       unsigned dimension,
                       ogdf::node source,
                       ogdf::node target,
                       bool directed = true);
    
private:
    const double epsilon_;
    
    struct Label {
        const Point * const point;
        ogdf::node n;
        const Label * const pred;
        bool mark_dominated;
        bool in_queue;
        
        inline Label(const Point *point, ogdf::node n, const Label *pred);
        inline Label(const Label &label);
        
        Label & operator=(const Label &label) = delete;
        
        ~Label() {
            delete point;
        }
    };
    
    struct LexLabelComp {
        bool operator()(const Label& l1, const Label& l2) {
            return LexPointComparator()(l2.point, l1.point);
        }
        
        bool operator()(const Label* l1, const Label* l2) {
            return LexPointComparator()(l2->point, l1->point);
        }
    };

};
    
EpSolverMartins::Label::
Label(const Point *point,
      ogdf::node n,
      const Label *pred)
:   point(point),
    n(n),
    pred(pred),
    mark_dominated(false),
    in_queue(false) {
}

EpSolverMartins::Label::
Label(const Label &label)
:   point(new Point(*label.point)),
    n(label.n),
    pred(label.pred),
    mark_dominated(label.mark_dominated),
    in_queue(label.in_queue) {
}

}

#endif /* MARTINS_H_ */
