#pragma once
/*

 * martins.h
 *
 *  Created on: 12.03.2013
 *      Author: fritz
 */

#ifndef WEIGHTED_MARTINS_B_H_
#define WEIGHTED_MARTINS_B_H_

#include <mco/basic/abstract_solver.h>

#include <setoper.h>
#include <cdd.h>

namespace mco {

class EpWeightedMartins : public AbstractSolver {

public:
	explicit EpWeightedMartins(double epsilon = 0)
    :   comp_leq_(epsilon, false) { }
    
	virtual void Solve(ogdf::Graph& graph,
                       std::function<const Point*(ogdf::edge)> weights,
                       unsigned dimension,
                       ogdf::node source,
                       ogdf::node target,
                       bool directed = true);
    
private:
    ComponentwisePointComparator comp_leq_;
    
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
    
    struct NodeEntry {
        std::list<Label*> label_set;
        std::list<unsigned> free_list;
        unsigned head;
        dd_MatrixPtr hull_matrix;
        
        inline void initialHullMatrix(unsigned dimension);
    };
    
    struct LexLabelComp {
        bool operator()(const Label& l1, const Label& l2) {
            return LexPointComparator()(l2.point, l1.point);
        }
        
        bool operator()(const Label* l1, const Label* l2) {
            return LexPointComparator()(l2->point, l1->point);
        }
    };

    bool ModifyLabelList(Label& new_label,
                         std::list<Label *>& label_set,
                         dd_MatrixPtr hull_matrix,
                         std::list<unsigned> free_list);

};
    
EpWeightedMartins::Label::
Label(const Point *point,
      ogdf::node n,
      const Label *pred)
:   point(point),
    n(n),
    pred(pred),
    mark_dominated(false),
    in_queue(false) {
}

EpWeightedMartins::Label::
Label(const Label &label)
:   point(new Point(*label.point)),
    n(label.n),
    pred(label.pred),
    mark_dominated(label.mark_dominated),
    in_queue(label.in_queue) {
}
    
void EpWeightedMartins::NodeEntry::
initialHullMatrix(unsigned dimension) {
    dd_MatrixPtr hull_matrix = dd_CreateMatrix(dimension + 1,
                                               dimension * 2);
    
    for(unsigned i = 0; i < dimension; ++i) {
        for(unsigned j = 0; j < dimension + 1; ++j) {
            dd_set_d(hull_matrix->matrix[i][j], i == j + 1 ? 1.0 : 0.0);
        }
        dd_set_d(hull_matrix->matrix[i][0], 0);
    }
    
}


}   // namespace mco

#endif /* WEIGHTED_MARTINS_H_ */
