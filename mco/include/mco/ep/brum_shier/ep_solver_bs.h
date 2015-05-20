#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_H_
#define BSSSA_H_

#include <mco/basic/abstract_solver.h>

namespace mco {

class EpSolverBS : public AbstractSolver<std::list<ogdf::edge>> {

    using heuristic_type = std::function<double(ogdf::node, unsigned)>;
    
public:
	EpSolverBS(double epsilon = 0)
    : epsilon_(epsilon) { }
    
	virtual void Solve(const ogdf::Graph& graph,
                       std::function<const Point*(const ogdf::edge)> costs,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed = true);

    void set_heuristic(heuristic_type heuristic) {
        heuristic_ = heuristic;
        use_heuristic_ = true;
    }

    void set_bounds(Point bounds) {
        bounds_ = std::move(bounds);
        use_bounds_ = true;
    }


private:
    const double epsilon_;
    unsigned dimension_;

    bool use_heuristic_ = false;
    heuristic_type heuristic_;
    bool use_bounds_ = false;
    Point bounds_;

    struct Label {
        const Point cost;
        const ogdf::node n;
        const ogdf::edge pred_edge;
        Label* pred_label;
        bool deleted = false;

        std::list<Label*> successors;

        Label(const Point cost,
              const ogdf::node n,
              const ogdf::edge p_edge,
              Label* p_label)
        :   cost(cost),
            n(n),
            pred_edge(p_edge),
            pred_label(p_label),
            deleted(false) { }

        Label(Label&& other)
        :   cost(std::move(other.cost)),
            n(other.n),
            pred_edge(other.pred_edge),
            pred_label(other.pred_label) {
        }
    };

    struct NodeEntry {
        using value_type = Label;

        bool in_queue;

        NodeEntry()
        :   in_queue(false),
            labels_(),
            labels_it_(labels_.begin()) { }

        NodeEntry(const NodeEntry& other)
        :   in_queue(false),
            labels_(),
            labels_it_(labels_.begin()) {
        }

        NodeEntry(NodeEntry&& other)
        :   in_queue(false),
            labels_(),
            labels_it_(labels_.begin()) {

        }

        void push_back(Label&& label) {
            if(labels_it_ == labels_.end()) {
                labels_it_ = labels_.insert(labels_it_, std::move(label));
            } else {
                labels_.push_back(std::move(label));
            }
        }

        std::list<Label>::iterator erase(std::list<Label>::iterator it) {
            if(it != labels_it_) {
                return labels_.erase(it);
            } else {
                return labels_it_ = labels_.erase(it);
            }
        }

        const std::list<Label>::iterator labels_it() {
            return labels_it_;
        }

        void proceed_labels_it() {
            labels_it_ = labels_.end();
        }

        bool has_new_labels() {
            return labels_it_ != labels_.end();
        }

        std::list<Label>& labels() {
            return labels_;
        }
    private:
        std::list<Label> labels_;
        std::list<Label>::iterator labels_it_;
    };

    bool check_domination(std::list<Label>& new_labels,
                          NodeEntry& neighbor_entry);

    void recursive_delete(Label& label);

    bool check_heuristic_prunable(const Label& label);


};

}

#endif /* BSSSA_H_ */
