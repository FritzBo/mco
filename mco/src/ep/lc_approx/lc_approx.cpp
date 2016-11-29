//
//  lc_approx.cpp
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#include <mco/ep/lc_approx/lc_approx.h>

using std::pair;
using std::list;
using std::vector;
using std::numeric_limits;

namespace mco {

bool LCApprox::
check_domination(vector<Label*>& new_labels,
                 NodeEntry& neighbor_entry)
{
    
    EqualityPointComparator eq;
    ComponentwisePointComparator leq(0, false);

    auto& neighbor_labels = neighbor_entry.labels();
    auto& neighbor_new_labels = neighbor_entry.new_labels();

    bool changed = false;

    bool dominated;

    for(auto new_label : new_labels)
    {

        dominated = false;

        unsigned check_label_it = 0;
        while(check_label_it < neighbor_labels.size()) {

            auto check_label = neighbor_labels[check_label_it];

            if(check_label->deleted) {

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

            } else if(eq(new_label->pos, check_label->pos) &&
                      new_label->sum < check_label->sum) {

                auto it = find(check_label->pred_label->succ_label.begin(),
                               check_label->pred_label->succ_label.end(),
                               check_label);

                assert(it != check_label->pred_label->succ_label.end());

                check_label->pred_label->succ_label.erase(it);

                recursive_delete(*check_label);

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

            } else if(leq(check_label->pos, new_label->pos)) {

                dominated = true;
                break;

            } else if(leq(new_label->pos, check_label->pos)) {

                auto it = find(check_label->pred_label->succ_label.begin(),
                               check_label->pred_label->succ_label.end(),
                               check_label);

                assert(it != check_label->pred_label->succ_label.end());

                check_label->pred_label->succ_label.erase(it);

                recursive_delete(*check_label);

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

            } else {
                ++check_label_it;
            }
        }

        if(!dominated) {

            check_label_it = 0;
            while(check_label_it < neighbor_new_labels.size()) {

                auto check_label = neighbor_new_labels[check_label_it];

                if(check_label->deleted) {

                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

                } else if(eq(new_label->pos, check_label->pos) &&
                          new_label->sum < check_label->sum) {

                    auto it = find(check_label->pred_label->succ_label.begin(),
                                   check_label->pred_label->succ_label.end(),
                                   check_label);

                    assert(it != check_label->pred_label->succ_label.end());

                    check_label->pred_label->succ_label.erase(it);

                    recursive_delete(*check_label);

                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

                } else if(leq(check_label->pos, new_label->pos)) {

                    dominated = true;
                    break;

                } else if(leq(new_label->pos, check_label->pos)) {

                    auto it = find(check_label->pred_label->succ_label.begin(),
                                   check_label->pred_label->succ_label.end(),
                                   check_label);

                    assert(it != check_label->pred_label->succ_label.end());
                    
                    check_label->pred_label->succ_label.erase(it);
                    
                    recursive_delete(*check_label);
                    
                    neighbor_entry.erase(neighbor_new_labels,
                                         check_label_it);

                } else {
                    ++check_label_it;
                }
            }
        }

        if(!dominated) {
            Label* pred = new_label->pred_label;
            neighbor_entry.push_back(new_label);
            pred->succ_label.push_back(new_label);

            changed = true;
        }
    }
    
    return changed;
}

void LCApprox::
recursive_delete(Label& label) {
    list<Label*> queue;
    queue.push_back(&label);

    while(!queue.empty()) {
        Label* curr = queue.front();
        queue.pop_front();

        curr->deleted = true;

        for(auto succ : curr->succ_label) {
            queue.push_back(succ);
        }
    }

}
    
void LCApprox::
Solve(const InstanceDescription& inst_desc) {

    const ForwardStar& graph =          *inst_desc.graph;
    cost_function_type cost_function =  inst_desc.cost_function;
    unsigned dimension =                inst_desc.dimension;
    const node source =                 inst_desc.source;
    const node target =                 inst_desc.target;
    const Point& epsilon =              *inst_desc.epsilon;

    epsilon_ = epsilon;
    dimension_ = dimension;
    min_e_ = Point(numeric_limits<double>::infinity(),
                   dimension_);
    
    for(auto e : graph.edges) {
        for(unsigned i = 0; i < dimension_; ++i) {
            min_e_[i] = std::min(min_e_[i], cost_function(e)[i]);
        }
    }

    list<Point> scaled_disj_bounds;
    Point scaled_bound(dimension_);

    FSNodeArray<NodeEntry> node_entries(graph);
    
    list<node> queue;
    
    {
        auto initial_label = new Label(Point(0.0, dimension),
                                       source,
                                       nulledge,
                                       nullptr,
                                       *this);
        
        node_entries[source].push_back(initial_label);
    }
    
    queue.push_back(source);
    node_entries[source].in_queue = true;
    
    while(!queue.empty()) {
        node current_node = queue.front();
        queue.pop_front();
        
        NodeEntry& current_node_entry = node_entries[current_node];
        current_node_entry.in_queue = false;
        
        if(current_node_entry.has_new_labels()) {

            auto& current_new_labels = current_node_entry.new_labels();
            
            for(auto current_edge : graph.adj_edges(current_node))
            {
                auto neighbor = graph.head(current_edge);
                
                auto& neighbor_entry = node_entries[neighbor];

                vector<Label*> new_labels;

                unsigned current_label_it = 0;
                while(current_label_it < current_node_entry.new_labels().size()) {

                    auto label = current_new_labels[current_label_it];

                    assert(label != nullptr);
                    assert(label->n == current_node);

                    if(!label->deleted) {

                        auto new_label = new Label(label->cost + cost_function(current_edge),
                                                   neighbor,
                                                   current_edge,
                                                   label,
                                                   *this);
                        

                        new_labels.push_back(new_label);
                    }

                    ++current_label_it;
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

    for(unsigned i = 0; i < node_entries[target].new_labels().size(); ++i) {
        auto label = node_entries[target].new_labels()[i];

        list<node> path;
        const Label* curr = label;
        if(!curr->deleted) {
            while(curr->n != source) {

                assert(graph.head(curr->pred_edge) == curr->n ||
                       graph.tail(curr->pred_edge) == curr->n);
                
                path.push_back(curr->n);
                curr = curr->pred_label;
                
            }
            
            assert(path.size() > 0);

            path.reverse();
            
            add_solution(path, label->cost);
        }
    }
    
}

}