/*
 * bsssa.cpp
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#include <mco/ep/brum_shier/ep_solver_bs.h>

#include <list>
#include <deque>
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

bool EpSolverBS::
check_domination(vector<Label*>& new_labels,
                 unsigned new_labels_end,
                 NodeEntry& neighbor_entry) {

    EqualityPointComparator eq;
    ComponentwisePointComparator leq(0, false);

    auto& neighbor_labels = neighbor_entry.labels();
    auto neighbor_labels_end = neighbor_entry.labels_end();

    auto& neighbor_new_labels = neighbor_entry.new_labels();
    auto neighbor_new_labels_end = neighbor_entry.new_labels_end();

    bool changed = false;

    bool dominated;

    unsigned new_label_it = 0;
    while(new_label_it < new_labels_end) {

        auto new_label = new_labels[new_label_it];

        assert(new_label != nullptr);

        dominated = false;

        unsigned check_label_it = 0;
        while(check_label_it < neighbor_labels_end) {

            auto check_label = neighbor_labels[check_label_it];

            if(check_label->deleted) {

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

#ifdef STATS
                ++deleted_labels_;
#endif

                --neighbor_labels_end;

            } else if(eq(new_label->cost, check_label->cost)) {

                remove_label(check_label);

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

#ifdef STATS
                ++deleted_labels_;
#endif

                --neighbor_labels_end;

#ifdef STATS
                label_compares_ += 1;
#endif

            } else if(leq(check_label->cost, new_label->cost)) {

                dominated = true;

#ifdef STATS
                label_compares_ += 2;
#endif

                break;

            } else if(leq(new_label->cost, check_label->cost)) {

                remove_label(check_label);

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

#ifdef STATS
                ++deleted_labels_;
#endif

                --neighbor_labels_end;

#ifdef STATS
                label_compares_ += 3;
#endif

            } else {
                ++check_label_it;

#ifdef STATS
                label_compares_ += 3;
#endif
            }
        }

        if(!dominated) {

            check_label_it = 0;
            while(check_label_it < neighbor_new_labels_end) {

                auto check_label = neighbor_new_labels[check_label_it];

                if(check_label->deleted) {

                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

#ifdef STATS
                    ++deleted_labels_;
#endif

                    --neighbor_new_labels_end;

                } else if(eq(new_label->cost, check_label->cost)) {

                    remove_label(check_label);

                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

#ifdef STATS
                    ++deleted_labels_;
#endif

                    --neighbor_new_labels_end;

#ifdef STATS
                    label_compares_ += 1;
#endif

                } else if(leq(check_label->cost, new_label->cost)) {

#ifdef STATS
                    label_compares_ += 2;
#endif

                    dominated = true;
                    break;

                } else if(leq(new_label->cost, check_label->cost)) {

                    remove_label(check_label);
                    
                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

#ifdef STATS
                    ++deleted_labels_;
#endif
                    
                    --neighbor_new_labels_end;

#ifdef STATS
                    label_compares_ += 3;
#endif

                } else {
                    ++check_label_it;

#ifdef STATS
                    label_compares_ += 3;
#endif

                }
            }
        }
        
        if(!dominated) {
            Label* pred = new_label->pred_label;
            neighbor_entry.push_back(new_label);
            ++neighbor_new_labels_end;
#if defined USE_TREE_DELETION || defined STATS
            pred->successors.push_back(new_label);
#endif

            changed = true;
        }
        else {
            delete new_label;
#ifdef STATS
            ++deleted_labels_;
#endif
        }
        
        ++new_label_it;
    }
    
    return changed;
}

#if defined USE_TREE_DELETION || defined STATS
void EpSolverBS::
recursive_delete(Label& label)
{
#if defined STATS && defined USE_TREE_DELETION
    recursive_deletions_ += 1;
#endif

    std::deque<Label*> queue;
    queue.push_back(&label);

    while(!queue.empty()) {
        Label* curr = queue.front();
        queue.pop_front();

#ifdef USE_TREE_DELETION
        curr->deleted = true;
#else
        curr->mark_recursive_deleted = true;
#endif

#if defined USE_TREE_DELETION && defined STATS
        ++deleted_tree_labels_;
#endif

        for(auto succ : curr->successors) {
            queue.push_back(succ);
        }
    }
}
#endif

bool EpSolverBS::check_heuristic_prunable(const Label& label) {

    for(unsigned i = 0; i < dimension_; ++i)
    {
        if(label.cost[i] + heuristic_(label.n, i) > bounds_[i])
        {
            return true;
        }
    }

    return false;
}


void EpSolverBS::Solve(const Graph& graph,
                       std::function<const Point*(const ogdf::edge)> weights,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed) {

    dimension_ = dimension;
    
    ring_buffer<node> queue(graph.numberOfNodes() + 1);
	NodeArray<NodeEntry> node_entries(graph);

    auto initial_label = new Label(Point(0.0, dimension),
                                   source,
                                   nullptr,
                                   nullptr);

    node_entries[source].push_back(initial_label);

	queue.push_back(source);
    node_entries[source].in_queue = true;

	while(!queue.empty()) {
		node current_node = queue.front();
        queue.pop_front();

//		cout << n << ": ";

		NodeEntry& current_node_entry = node_entries[current_node];
        current_node_entry.in_queue = false;

        if(current_node_entry.has_new_labels()) {

            auto& current_new_labels = current_node_entry.new_labels();

            for(auto adj : current_node->adjEntries) {
                
                edge current_edge = adj->theEdge();
                
                if(current_edge->isSelfLoop())
                {
                    continue;
                }

                if(directed && current_edge->target() == current_node)
                {
                    continue;
                }

#ifdef STATS
                arc_pushes_ += 1;
#endif

                node neighbor = current_edge->opposite(current_node);
                auto& neighbor_entry = node_entries[neighbor];

    //			cout << neighbor << ", ";

                vector<Label*> new_labels(current_node_entry.new_labels_end(), nullptr);

                unsigned size = 0;

                unsigned current_label_it = 0;
                while(current_label_it < current_node_entry. new_labels_end()) {

                    auto label = current_new_labels[current_label_it];

#ifdef STATS
                    if(label->mark_recursive_deleted)
                    {
                        ++touched_recursively_deleted_label_;
                    }
#endif

                    assert(label != nullptr);
                    assert(label->n == current_node);

                    if(!label->deleted) {

                        auto new_label = new Label(label->cost + *weights(current_edge),
                                                   neighbor,
                                                   current_edge,
                                                   label);

                        if(!use_bounds_ || !use_heuristic_ ||
                           !check_heuristic_prunable(*new_label)) {

                            new_labels[size] = new_label;
                            ++size;
                        } else {
                            delete new_label;
#ifdef STATS
                            ++deleted_labels_;
#endif
                        }
                    }

                    ++current_label_it;
                }

                if(!new_labels.empty()) {


                    bool changed = check_domination(new_labels,
                                                    size,
                                                    neighbor_entry);

                    if(!node_entries[neighbor].in_queue &&
                       changed &&
                       neighbor != target &&
                       neighbor != source) {

                        queue.push_back(neighbor);
                        node_entries[neighbor].in_queue = true;
                    }

                }
            }

            current_node_entry.proceed_labels_it();

        }

    }

    auto& target_entry = node_entries[target];
    for(unsigned i = 0; i < target_entry.new_labels_end(); ++i) {
        auto label = target_entry.new_labels()[i];
        list<edge> path;
        auto curr = label;
        if(!curr->deleted) {
            while(curr->n != source) {

                assert(curr->pred_edge->source() == curr->n ||
                       curr->pred_edge->target() == curr->n);

                path.push_back(curr->pred_edge);
                curr = curr->pred_label;

            }

            assert(path.size() > 0);

            path.reverse();

            add_solution(path, label->cost);
        }
    }
}

}
