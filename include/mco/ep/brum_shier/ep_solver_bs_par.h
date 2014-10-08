#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_PAR_H_
#define BSSSA__PAR_H_

#include <mco/basic/abstract_solver.h>

#include <thread>

namespace mco {

class EpSolverBSPar : public AbstractSolver<std::list<ogdf::edge>> {

    using heuristic_type = std::function<double(ogdf::node, unsigned)>;
    using cost_type = std::function<const Point*(const ogdf::edge)>;
    
public:

	virtual void Solve(const ogdf::Graph& graph,
                       cost_type costs,
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

    void set_no_threads(unsigned no) {
        no_threads_ = no;
    }


private:
    unsigned dimension_;
    bool directed_;
    ogdf::node source_;
    ogdf::node target_;
    cost_type weights_;

    unsigned no_threads_ = 1;
    std::atomic<unsigned> pushing_threads_;

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

    ogdf::NodeArray<NodeEntry>* node_entries_;

    std::list<ogdf::node>* queue_;
    std::mutex queue_mutex_;

    std::mutex neighborhood_lock_;
    ogdf::NodeArray<bool>* neighborhood_locks_;
    bool lock_neighborhood(ogdf::node n);
    void unlock_neighborhood(ogdf::node n);

    void thread_worker();

    bool check_domination(std::list<Label>& new_labels,
                          NodeEntry& neighbor_entry);

    bool check_heuristic_prunable(const Label& label);

    void push_labels(ogdf::node current_node);


};

}

#endif /* BSSSA_H_ */
