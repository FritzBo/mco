/*
 * bsssa.cpp
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#include <mco/ep/brum_shier/ep_solver_bs.h>

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

            if(check_label.deleted) {

                check_label_it = neighbor_entry.erase(check_label_it);

            } else if(leq(check_label.cost, new_label.cost)) {

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

bool EpSolverBS::check_heuristic_prunable(const Label& label) {

    Point heuristic_cost(dimension_);
    for(unsigned i = 0; i < dimension_; ++i) {
        heuristic_cost[i] = label.cost[i] + heuristic_(label.n, i);
    }

    return !ComponentwisePointComparator(0, false)(heuristic_cost, bounds_);
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
        Label initial_label(Point(0.0, dimension),
                            source,
                            nullptr,
                            nullptr);

        node_entries[source].push_back(std::move(initial_label));
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

                list<Label> new_labels;

                for(auto current_label_it = current_node_entry.labels_it();
                    current_label_it != current_node_entry.labels().end();
                    ++current_label_it) {

                    Label& label = *current_label_it;

                    assert(label.n == current_node);

                    if(!label.deleted) {

                        Label new_label(label.cost + *weights(current_edge),
                                        neighbor,
                                        current_edge,
                                        &label);

                        if(!use_bounds_ || !use_heuristic_ ||
                           !check_heuristic_prunable(new_label)) {

                            new_labels.push_back(std::move(new_label));
                        }
                    }
                }

                if(!new_labels.empty()) {


                    bool changed = check_domination(new_labels,
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
    
    for(auto& label : node_entries[target].labels()) {
        list<edge> path;
        const Label* curr = &label;
        if(!curr->deleted) {
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
    }
}

}
