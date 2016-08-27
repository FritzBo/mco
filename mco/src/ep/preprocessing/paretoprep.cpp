//
//  paretoprep.cpp
//  mco
//
//  Created by Fritz Bökler on 18.02.16.
//
//

#include <mco/ep/preprocessing/paretoprep.h>

using std::numeric_limits;
using std::make_pair;
using std::vector;

#include <ogdf/basic/PriorityQueue.h>

using ogdf::PrioritizedMapQueue;

namespace mco {

void ParetoPrep::preprocess(const ForwardStar& reverse_star,
                            const ReverseEdgeArray<Point>& weights,
                            unsigned no_objectives,
                            node source,
                            node target)
{

    ComponentwisePointComparator comp;

    std::vector<Point> upper_bounds;
    FSNodeArray<NodeEntry> node_entries(reverse_star, no_objectives);

    for(node n : reverse_star.nodes)
    {
        for(unsigned i = 0; i < no_objectives; ++i)
        {
            node_entries[n].temp_lower_bounds[i] = numeric_limits<double>::infinity();
            node_entries[n].predecessor_edges[i] = nulledge;
        }
    }

    for(unsigned i = 0; i < no_objectives; ++i)
    {
        node_entries[target].temp_lower_bounds[i] = 0;
        node_entries[target].predecessor_edges[i] = target;
    }

    PrioritizedMapQueue<node, double> queue;
    queue.push(target, 0);

    while(!queue.empty())
    {
        node current_node = queue.topElement();
        queue.pop();

        NodeEntry& current_entry = node_entries[current_node];

        bool temp_dominated_node = false;
        for(auto& ub : upper_bounds)
        {
            if(comp(ub, current_entry.temp_lower_bounds))
            {
                temp_dominated_node = true;
            }
        }

        if(temp_dominated_node)
        {
            continue;
        }

        for(auto edge : reverse_star.adj_edges(current_node))
        {
            node in_neighbor = reverse_star.head(edge);
            NodeEntry& neighbor_entry = node_entries[in_neighbor];

            bool changed = false;

            for(unsigned i = 0; i < no_objectives; ++i)
            {
                double neighbor_lb = neighbor_entry.temp_lower_bounds[i];
                double edge_cost = weights[edge][i];
                double own_lb = current_entry.temp_lower_bounds[i];

                if(own_lb + edge_cost < neighbor_lb)
                {
                    neighbor_entry.temp_lower_bounds[i] = own_lb + edge_cost;
                    neighbor_entry.predecessor_edges[i] = edge;

                    changed = true;

                    if(in_neighbor == source)
                    {
                        Point new_upper_bound(no_objectives);
                        construct_path(reverse_star,
                                       weights,
                                       node_entries,
                                       source,
                                       target,
                                       i,
                                       new_upper_bound);

                        upper_bounds.push_back(std::move(new_upper_bound));
                    }

                }
            }

            if(changed)
            {
                double sum = 0;
                for(unsigned i = 0; i < no_objectives; ++i)
                {
                    sum += neighbor_entry.temp_lower_bounds[i];
                }

                if(queue.contains(in_neighbor))
                {
                    queue.decrease(in_neighbor, sum);
                }
                else
                {
                    queue.push(in_neighbor, sum);
                }
            }
        }
    }

    ideal_point_ = std::move(node_entries[source].temp_lower_bounds);

#ifndef NDEBUG
    std::cout << "ideal point:" << std::endl;
    std::cout << node_entries[source].temp_lower_bounds << std::endl;

    std::cout << "paths:" << std::endl;
    for(auto& p : upper_bounds)
    {
        std::cout << p << std::endl;
    }

    unsigned deleted_nodes = 0;
    for(node n : reverse_star.nodes)
    {
        bool not_deleted = true;
        for(unsigned i = 0; i < no_objectives; ++i)
        {
            if(node_entries[n].predecessor_edges[i] == nulledge)
            {
                not_deleted = false;
            }
        }

        if(not_deleted == false)
        {
            deleted_nodes += 1;
        }
    }
    std::cout << "deleted nodes: " << deleted_nodes << std::endl;
#endif
}

void ParetoPrep::construct_path(const ForwardStar& reverse_star,
                                const ReverseEdgeArray<Point>& weights,
                                const FSNodeArray<NodeEntry>& node_entries,
                                node source,
                                node target,
                                unsigned objective,
                                Point& new_upper_bound)
{
    node current_node = source;

    while(current_node != target)
    {
        auto& current_entry = node_entries[current_node];
        edge e = current_entry.predecessor_edges[objective];

        new_upper_bound += weights[e];

        current_node = reverse_star.tail(e);
    }
}

}