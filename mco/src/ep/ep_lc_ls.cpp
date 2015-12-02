/*
 * bsssa.cpp
 *
 *  Created on: 1.12.2015
 *      Author: fritz
 */

#include <mco/ep/ep_lc_ls.h>

#include <deque>
#include <list>
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include <functional>

using std::deque;
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

bool EpLcLs::
check_domination(Label* new_label,
                 NodeEntry& neighbor_entry)
{

    EqualityPointComparator eq;
    ComponentwisePointComparator leq(0, false);

    auto& neighbor_labels = neighbor_entry.labels();
    auto neighbor_labels_end = neighbor_entry.labels_end();

    bool changed = false;

    bool dominated = false;

    unsigned check_label_it = 0;
    while(check_label_it < neighbor_labels_end)
    {

        auto check_label = neighbor_labels[check_label_it];

        if(check_label->deleted)
        {

            neighbor_entry.erase(check_label_it);

            --neighbor_labels_end;

        }
        else if(eq(new_label->cost, check_label->cost))
        {

            auto it = find(check_label->pred_label->successors.begin(),
                           check_label->pred_label->successors.end(),
                           check_label);

            assert(it != check_label->pred_label->successors.end());

            check_label->pred_label->successors.erase(it);

            recursive_delete(*check_label);

            neighbor_entry.erase(check_label_it);

            --neighbor_labels_end;

        }
        else if(leq(check_label->cost, new_label->cost))
        {

            dominated = true;
            break;

        }
        else if(leq(new_label->cost, check_label->cost))
        {

            auto it = find(check_label->pred_label->successors.begin(),
                           check_label->pred_label->successors.end(),
                           check_label);

            assert(it != check_label->pred_label->successors.end());

            check_label->pred_label->successors.erase(it);

            recursive_delete(*check_label);

            neighbor_entry.erase(check_label_it);

            --neighbor_labels_end;

        }
        else
        {
            ++check_label_it;
        }
    }

    if(!dominated)
    {
        Label* pred = new_label->pred_label;
        neighbor_entry.push_back(new_label);
        pred->successors.push_back(new_label);

        new_label->in_node_list = true;
        
        changed = true;
    }

    return changed;
}

void EpLcLs::recursive_delete(Label& label)
{
    // deque superior
    deque<Label*> queue;
    queue.push_back(&label);

    while(!queue.empty())
    {
        Label* curr = queue.front();
        queue.pop_front();

        curr->deleted = true;

        for(auto succ : curr->successors)
        {
            queue.push_back(succ);
        }
    }
}

bool EpLcLs::check_heuristic_prunable(const Label& label)
{

    Point heuristic_cost(dimension_);
    for(unsigned i = 0; i < dimension_; ++i)
    {
        heuristic_cost[i] = label.cost[i] + heuristic_(label.n, i);
    }

    return !ComponentwisePointComparator(0, false)(heuristic_cost, bounds_);
}


void EpLcLs::Solve(const Graph& graph,
                   std::function<const Point*(const ogdf::edge)> weights,
                   unsigned dimension,
                   const ogdf::node source,
                   const ogdf::node target,
                   bool directed)
{

    dimension_ = dimension;

    /*  deque is the superior datastructure here,
        because the queue might grow very large */

	deque<Label*> queue;
	NodeArray<NodeEntry> node_entries(graph);

    auto initial_label = new Label(Point(0.0, dimension),
                                   source,
                                   nullptr,
                                   nullptr);

    node_entries[source].push_back(initial_label);
    queue.push_back(initial_label);
    initial_label->in_queue = true;

	while(!queue.empty())
    {
		auto current_label = queue.front();
        queue.pop_front();
        current_label->in_queue = false;

        if(!current_label->deleted)
        {

            node current_node = current_label->n;

            if(current_node == target)
            {
                continue;
            }

            for(auto adj : current_node->adjEdges)
            {

                edge current_edge = adj->theEdge();
                    
                if(current_edge->isSelfLoop())
                {
                    continue;
                }

                if(directed && current_edge->target() == current_node)
                {
                    continue;
                }

                node neighbor = current_edge->opposite(current_node);
                auto& neighbor_entry = node_entries[neighbor];

                if(neighbor == source) {
                    continue;
                }

                auto new_label = new Label(current_label->cost + *weights(current_edge),
                                           neighbor,
                                           current_edge,
                                           current_label);

                if(!use_bounds_ || !use_heuristic_ ||
                   !check_heuristic_prunable(*new_label))
                {

                    bool changed = check_domination(new_label,
                                                    neighbor_entry);

                    if(changed)
                    {

                        queue.push_back(new_label);
                        new_label->in_queue = true;
                    }
                    else
                    {
                        delete new_label;
                    }

                }

            }

        }
        else if(!current_label->in_node_list)
        {
            delete current_label;
        }

    }
    
    for(auto label : node_entries[target].labels())
    {
        list<edge> path;
        auto curr = label;
        if(!curr->deleted)
        {
            while(curr->n != source)
            {
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
