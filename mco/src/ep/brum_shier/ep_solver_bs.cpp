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

                --neighbor_labels_end;

            } else if(eq(new_label->cost, check_label->cost)) {

//                check_label->pred_label->erase_successor(check_label);

                auto it = find(check_label->pred_label->successors.begin(),
                               check_label->pred_label->successors.end(),
                               check_label);

                assert(it != check_label->pred_label->successors.end());

                check_label->pred_label->successors.erase(it);

                recursive_delete(*check_label);

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

                --neighbor_labels_end;

            } else if(leq(check_label->cost, new_label->cost)) {

                dominated = true;
                break;

            } else if(leq(new_label->cost, check_label->cost)) {

//                check_label->pred_label->erase_successor(check_label);

                auto it = find(check_label->pred_label->successors.begin(),
                               check_label->pred_label->successors.end(),
                               check_label);

                assert(it != check_label->pred_label->successors.end());

                check_label->pred_label->successors.erase(it);


                recursive_delete(*check_label);

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

                --neighbor_labels_end;

            } else {
                ++check_label_it;
            }
        }

        if(!dominated) {

            check_label_it = 0;
            while(check_label_it < neighbor_new_labels_end) {

                auto check_label = neighbor_new_labels[check_label_it];

                if(check_label->deleted) {

                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

                    --neighbor_new_labels_end;

                } else if(eq(new_label->cost, check_label->cost)) {

//                    check_label->pred_label->erase_successor(check_label);

                    auto it = find(check_label->pred_label->successors.begin(),
                                   check_label->pred_label->successors.end(),
                                   check_label);

                    assert(it != check_label->pred_label->successors.end());

                    check_label->pred_label->successors.erase(it);

                    recursive_delete(*check_label);

                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

                    --neighbor_new_labels_end;

                } else if(leq(check_label->cost, new_label->cost)) {

                    dominated = true;
                    break;

                } else if(leq(new_label->cost, check_label->cost)) {
                    
//                    check_label->pred_label->erase_successor(check_label);

                    auto it = find(check_label->pred_label->successors.begin(),
                                   check_label->pred_label->successors.end(),
                                   check_label);

                    assert(it != check_label->pred_label->successors.end());

                    check_label->pred_label->successors.erase(it);

                    
                    recursive_delete(*check_label);
                    
                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);
                    
                    --neighbor_new_labels_end;
                    
                } else {
                    ++check_label_it;
                }
            }
        }
        
        if(!dominated) {
            Label* pred = new_label->pred_label;
            neighbor_entry.push_back(new_label);
            ++neighbor_new_labels_end;
            pred->successors.push_back(new_label);
//            pred->push_successor(new_label);

            changed = true;
        }
        else {
            delete new_label;
        }
        
        ++new_label_it;
    }
    
    return changed;
}

void EpSolverBS::
recursive_delete(Label& label)
{
    std::deque<Label*> queue;
    queue.push_back(&label);

    while(!queue.empty()) {
        Label* curr = queue.front();
        queue.pop_front();

        curr->deleted = true;

        for(auto succ : curr->successors) {
            queue.push_back(succ);
        }
    }
}

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
    
	list<node> queue;
	NodeArray<NodeEntry> node_entries(graph);

    {
        auto initial_label = new Label(Point(0.0, dimension),
                                       source,
                                       nullptr,
                                       nullptr);

        node_entries[source].push_back(initial_label);
    }

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

            for(auto adj : current_node->adjEdges) {
                
                edge current_edge = adj->theEdge();
                
                if(current_edge->isSelfLoop()) {
                    continue;
                }

                if(directed && current_edge->target() == current_node) {
                    continue;
                }

                node neighbor = current_edge->opposite(current_node);
                auto& neighbor_entry = node_entries[neighbor];

    //			cout << neighbor << ", ";

                vector<Label*> new_labels(current_node_entry.new_labels_end(), nullptr);

                unsigned size = 0;

                unsigned current_label_it = 0;
                while(current_label_it < current_node_entry. new_labels_end()) {

                    auto label = current_new_labels[current_label_it];

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
    
    for(auto label : node_entries[target].new_labels()) {
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
