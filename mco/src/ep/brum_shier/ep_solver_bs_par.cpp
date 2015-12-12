/*
 * bsssa.cpp
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#include <mco/ep/brum_shier/ep_solver_bs_par.h>

#include <list>
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include <functional>

using std::list;
using std::vector;
using std::cout;
using std::endl;
using std::function;
using std::pair;
using std::thread;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::node;
using ogdf::Graph;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using ogdf::AdjElement;

#include <mco/basic/point.h>
#include <mco/ep/basic/ep_instance.h>
#include <mco/basic/utility.h>

namespace mco {

bool EpSolverBSPar::
check_domination(list<Label>& new_labels,
                 NodeEntry& neighbor_entry) {

    EqualityPointComparator eq;
    ComponentwisePointComparator leq(0, false);

    bool changed = false;

    list<Label>& neighbor_labels = neighbor_entry.labels();

    bool dominated;

    for(auto& new_label : new_labels) {
        dominated = false;

        auto check_label_it = neighbor_labels.begin();
        while(check_label_it != neighbor_labels.end()) {

            Label& check_label = *check_label_it;

            if(leq(check_label.cost, new_label.cost)) {

                dominated = true;
                break;

            } else if(leq(new_label.cost, check_label.cost)) {

                check_label_it = neighbor_entry.erase(check_label_it);

                changed = true;

            } else {
                ++check_label_it;
            }
        }
        
        if(!dominated) {
            neighbor_entry.push_back(std::move(new_label));
            changed = true;
        }
        
    }
    
    return changed;
}

bool EpSolverBSPar::check_heuristic_prunable(const Label& label) {

    Point heuristic_cost(dimension_);
    for(unsigned i = 0; i < dimension_; ++i) {
        heuristic_cost[i] = label.cost[i] + heuristic_(label.n, i);
    }

    return !ComponentwisePointComparator(0, false)(heuristic_cost, bounds_);
}


void EpSolverBSPar::Solve(const Graph& graph,
                       cost_type weights,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed) {

    dimension_ = dimension;
    directed_ = directed;
    source_ = source;
    target_ = target;
    weights_ = weights;
    
	queue_ = new list<node>;
	node_entries_ = new NodeArray<NodeEntry>(graph);

    neighborhood_locks_ = new NodeArray<bool>(graph, false);

    list<node>& queue = *queue_;
    NodeArray<NodeEntry>& node_entries = *node_entries_;

    {
        Label initial_label(Point(0.0, dimension),
                            source,
                            nullptr,
                            nullptr);

        node_entries[source].push_back(std::move(initial_label));
    }

	queue.push_back(source);
    node_entries[source].in_queue = true;
    pushing_threads_.store(0);

    list<thread*> thread_list;
    for(unsigned i = 0; i < no_threads_; ++i) {
        thread* new_thread;
        new_thread = new thread(&EpSolverBSPar::thread_worker, this);
        thread_list.push_back(new_thread);
    }

    while(!thread_list.empty()) {
        thread_list.front()->join();
        thread_list.pop_front();
    }


    for(auto& label : node_entries[target].labels()) {
        list<edge> path;
        const Label* curr = &label;
        while(curr->n != source) {

            assert(curr->pred_edge->source() == curr->n ||
                   curr->pred_edge->target() == curr->n);

            path.push_back(curr->pred_edge);
            curr = curr->pred_label;

        }

        assert(path.size() > 0);

        path.reverse();

        add_solution(path, label.cost);
    }

    delete node_entries_;
    delete queue_;
    delete neighborhood_locks_;
}

void EpSolverBSPar::thread_worker() {
    list<node>& queue = *queue_;

    queue_mutex_.lock();
    while(!queue.empty() || pushing_threads_.load() > 0) {

        if(!queue.empty()) {

            bool still_locked = true;
            auto queue_it = queue.begin();
            while(queue_it != queue.end()) {
                auto current_node = *queue_it;

                if(lock_neighborhood(current_node)) {
                    queue.erase(queue_it);
                    queue_mutex_.unlock();
                    still_locked = false;

                    pushing_threads_++;
                    push_labels(current_node);
                    pushing_threads_--;
                    unlock_neighborhood(current_node);
                    
                    break;
                }

                queue_it++;
            }

            if(still_locked) {
                queue_mutex_.unlock();
            }

        } else {
            queue_mutex_.unlock();
        }

        while(queue.empty() && pushing_threads_.load() > 0) {
            std::this_thread::yield();
        }

        queue_mutex_.lock();
    }

    queue_mutex_.unlock();

}

void EpSolverBSPar::push_labels(node current_node) {
    list<node>& queue = *queue_;
    NodeArray<NodeEntry>& node_entries = *node_entries_;

    NodeEntry& current_node_entry = node_entries[current_node];
    current_node_entry.in_queue = false;

    if(current_node_entry.has_new_labels()) {


        for(auto adj : current_node->adjEntries) {

            edge current_edge = adj->theEdge();

            if(current_edge->isSelfLoop()) {
                continue;
            }

            if(directed_ && current_edge->target() == current_node) {
                continue;
            }

            node neighbor = current_edge->opposite(current_node);
            auto& neighbor_entry = node_entries[neighbor];

            //			cout << neighbor << ", ";

            list<Label> new_labels;

            for(auto current_label_it = current_node_entry.labels_it();
                current_label_it != current_node_entry.labels().end();
                ++current_label_it) {

                Label& label = *current_label_it;

                assert(label.n == current_node);

                Label new_label(label.cost + *weights_(current_edge),
                                neighbor,
                                current_edge,
                                &label);

                if(!use_bounds_ || !use_heuristic_ ||
                   !check_heuristic_prunable(new_label)) {

                    new_labels.push_back(std::move(new_label));
                }
            }

            if(!new_labels.empty()) {


                bool changed = check_domination(new_labels,
                                                neighbor_entry);

                if(!node_entries[neighbor].in_queue &&
                   changed &&
                   neighbor != target_ &&
                   neighbor != source_) {

                    queue_mutex_.lock();
                    queue.push_back(neighbor);
                    queue_mutex_.unlock();
                    
                    node_entries[neighbor].in_queue = true;
                }

            }
        }

        current_node_entry.proceed_labels_it();

    }

}

bool EpSolverBSPar::lock_neighborhood(node n) {
    neighborhood_lock_.lock();

    auto& neighborhood_locks = *neighborhood_locks_;

    if(neighborhood_locks[n]) {
        neighborhood_lock_.unlock();
        return false;
    }

    for(auto adj : n->adjEntries) {
        auto edge = adj->theEdge();
        auto neighbor = edge->opposite(n);

        if(neighborhood_locks[neighbor]) {
            neighborhood_lock_.unlock();
            return false;
        }

    }

    neighborhood_locks[n] = true;

    for(auto adj : n->adjEntries) {
        auto edge = adj->theEdge();
        auto neighbor = edge->opposite(n);

        neighborhood_locks[neighbor] = true;
        
    }

    neighborhood_lock_.unlock();

    return true;
}

void EpSolverBSPar::unlock_neighborhood(node n) {
    neighborhood_lock_.lock();
    auto& neighborhood_locks = *neighborhood_locks_;

    neighborhood_locks[n] = false;

    for(auto adj : n->adjEntries) {
        auto edge = adj->theEdge();
        auto neighbor = edge->opposite(n);

        neighborhood_locks[neighbor] = false;

    }

    neighborhood_lock_.unlock();
}

}
