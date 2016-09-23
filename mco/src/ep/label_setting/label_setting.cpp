//
//  label_selection.cpp
//  mco
//
//  Created by Fritz Bökler on 16.09.16.
//
//

#include <mco/ep/label_setting/label_setting.h>

#include <queue>

using std::function;
using std::vector;
using std::priority_queue;

using mco::MaxAryHeap;

namespace mco {

bool MOSPLabelSetting::non_dominated(Point new_cost,
                                     NodeEntry& target_node_entry)
{
    ComponentwisePointComparator comp(0, false);

    // Check if a closed label dominates the new cost
    for(auto& target_label : target_node_entry.closed_labels())
    {
        if(comp(target_label.cost, new_cost))
        {
            return false;
        }
    }

    // Check if an open label dominates the new cost
    // or if the new cost dominates one of the open labels
    bool dominated_another = false;
    for(size_t i = 0; i < target_node_entry.open_labels().size(); ++i)
    {
        auto& target_label = target_node_entry.open_labels()[i];

        if(target_label.deleted)
        {
            continue;
        }
        else if(!dominated_another && comp(target_label.cost, new_cost))
        {
            return false;
        }
        else if(comp(new_cost, target_label.cost))
        {
            target_label.deleted = true;
            dominated_another = true;
        }
    }

    return true;
}

void MOSPLabelSetting::Solve(InputDefinition& def)
{
    const ForwardStar& graph = *def.graph;
    function<const Point*(const edge)> costs = def.costs;
    unsigned dimension = def.dimension;
    const node source = def.source;
    const node target = def.target;

    FSNodeArray<NodeEntry> node_entries(graph);

    vector<node> queue;
    NodeComparator cmp(node_entries);

    {
        Label initial_label(Point(dimension),
                            source,
                            nulledge,
                            0,
                            node_entries[source]);

        node_entries[source].push_open(std::move(initial_label));
        node_entries[source].in_queue = true;

        queue.push_back(source);
        push_ary_heap<D_nodes>(queue.begin(),
                               queue.end(),
                               cmp);
    }

    Point new_cost(dimension);

    while(!queue.empty())
    {
        auto current_node = queue.front();
        pop_ary_heap<D_nodes>(queue.begin(),
                              queue.end(),
                              cmp);
        queue.pop_back();


        auto& current_node_entry = node_entries[current_node];

        auto current_label = std::move(current_node_entry.top_open());
        current_node_entry.pop_open();

#ifndef NDEBUG
        std::cout << current_node << " (" << current_label.sum << ")" << std::endl;
#endif

        if(current_node_entry.has_open())
        {
            queue.push_back(current_node);
            push_ary_heap<D_nodes>(queue.begin(),
                                   queue.end(),
                                   cmp);
        } else {
            current_node_entry.in_queue = false;
        }

        if(current_label.deleted)
        {
            continue;
        }

        if(current_node == target)
        {
            current_node_entry.close_label(std::move(current_label));
            continue;
        }

        for(auto e : graph.adj_edges(current_node))
        {
            auto neighbor = graph.head(e);
            auto& neighbor_node_entry = node_entries[neighbor];

#ifndef NDEBUG
            std::cout << "-> " << neighbor;
#endif

            new_cost = current_label.cost + *costs(e);

            if(non_dominated(new_cost, neighbor_node_entry))
            {
                Label new_label(std::move(new_cost),
                                neighbor,
                                e,
                                current_label.label_id,
                                neighbor_node_entry);

#ifndef NDEBUG
                std::cout << " (" << new_label.sum << ")" << std::endl;
#endif

                neighbor_node_entry.push_open(std::move(new_label));

                if(neighbor != source &&
                   !neighbor_node_entry.in_queue)
                {
                    queue.push_back(neighbor);
                    push_ary_heap<D_nodes>(queue.begin(),
                                           queue.end(),
                                           cmp);
                    neighbor_node_entry.in_queue = true;
                }

            }
#ifndef NDEBUG
            else
            {
                std::cout << std::endl;
            }
#endif
        }

        current_node_entry.close_label(std::move(current_label));
    }
    
    std::cout << node_entries[target].closed_labels().size() << std::endl;

//    for(auto& l : node_entries[target].closed_labels())
//    {
//        std::cout << l.cost << std::endl;
//    }
}

}