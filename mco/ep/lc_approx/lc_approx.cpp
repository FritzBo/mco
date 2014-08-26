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

using ogdf::node;
using ogdf::edge;
using ogdf::Graph;
using ogdf::NodeArray;

namespace mco {

bool LCApprox::
check_domination(list<Label>& new_labels,
                 NodeEntry& neighbor_entry) {
    
    EqualityPointComparator eq;
    ComponentwisePointComparator leq(0, false);
    
    bool changed = false;
    
    list<Label>& neighbor_labels = neighbor_entry.labels();
    
    bool dominated;
    
    for(auto new_label : new_labels) {
        dominated = false;
        
        auto check_label_it = neighbor_labels.begin();
        while(check_label_it != neighbor_labels.end()) {
            
            const Label& check_label = *check_label_it;
            
            if(leq(check_label.pos, new_label.pos)) {
                dominated = true;
                break;
            }
            
            if(leq(new_label.pos, check_label.pos)) {
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
    
bool LCApprox::check_heuristic_prunable(const Label& label,
                                        const Point& bounds) {
    
    Point heuristic_cost(dimension_);
    for(unsigned i = 0; i < dimension_; ++i) {
        heuristic_cost[i] = label.cost[i] + heuristic_(label.n, i);
    }
    
    compute_pos(heuristic_cost, heuristic_cost);
    
    return !ComponentwisePointComparator(0, false)(heuristic_cost, bounds);
}
    
void LCApprox::
Solve(const Graph& graph,
      cost_function_type cost_function,
      unsigned dimension,
      const node source,
      const node target,
      bool directed,
      const Point& epsilon) {
    
    epsilon_ = epsilon;
    dimension_ = dimension;
    min_e_ = Point(numeric_limits<double>::infinity(),
                   dimension_);
    
    for(auto e : graph.edges) {
        for(unsigned i = 0; i < dimension_; ++i) {
            min_e_[i] = min(min_e_[i], cost_function(e)[i]);
        }
    }
    
    Point scaled_bounds(dimension_);
    
    if(use_bounds_) {
        compute_pos(bounds_, scaled_bounds);
    }
    
    NodeArray<NodeEntry> node_entries(graph);
    
    list<node> queue;
    
    {
        Label initial_label(Point(0.0, dimension),
                            source,
                            nullptr,
                            nullptr,
                            *this);
        
        node_entries[source].push_back(std::move(initial_label));
    }
    
    queue.push_back(source);
    node_entries[source].in_queue = true;
    
    while(!queue.empty()) {
        node current_node = queue.front();
        queue.pop_front();
        
        NodeEntry& current_node_entry = node_entries[current_node];
        current_node_entry.in_queue = false;
        
        if(current_node_entry.has_new_labels()) {
            
            for(auto adj : current_node->adjEdges) {
                edge current_edge = adj->theEdge();
                
                if(directed && current_edge->target() == current_node) {
                    continue;
                }
                
                auto neighbor = current_edge->target() != current_node ? current_edge->target() : current_edge->source();
                
                auto& neighbor_entry = node_entries[neighbor];
                
                list<Label> new_labels;
                
                for(auto current_label_it = current_node_entry.labels_it();
                    current_label_it != current_node_entry.labels().end();
                    ++current_label_it) {
                    
                    Label new_label(current_label_it->cost + cost_function(current_edge),
                                    neighbor,
                                    current_edge,
                                    &*current_label_it,
                                    *this);
                    
                    if(!use_bounds_ || !use_heuristic_ ||
                       !check_heuristic_prunable(new_label, scaled_bounds)) {
                        
                        new_labels.push_back(std::move(new_label));
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
                    }
                    
                }
                
                
            }

            current_node_entry.proceed_labels_it();
        }
        
    }
    
	for(auto label : node_entries[target].labels()) {
        list<edge> path;
        const Label* curr = &label;
        while(curr->n != source) {
            
            path.push_back(curr->pred_edge);
            curr = curr->pred_label;
            
        }
        
        assert(path.size() > 0);
        
        path.reverse();
        
        add_solution(path, label.cost);
    }
    
}

}