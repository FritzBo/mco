//
//  lc_approx.cpp
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#include <mco/ep/lc_approx/lc_approx.h>

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
    
    const list<Label>& neighbor_labels = neighbor_entry.labels();
    
    bool dominated;
    
    for(auto new_label : new_labels) {
        dominated = false;
        
        auto check_label_it = neighbor_labels.cbegin();
        while(check_label_it != neighbor_labels.cend()) {
            
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
    
void LCApprox::
Solve(const Graph& graph,
      cost_function_type cost_function,
      unsigned dimension,
      const node source,
      const node target,
      bool directed,
      double epsilon) {
    
    epsilon_ = epsilon;
    dimension_ = dimension;
    min_e_ = Point(numeric_limits<double>::infinity(),
                   dimension_);
    
    for(auto e : graph.edges) {
        for(unsigned i = 0; i < dimension_; ++i) {
            min_e_[i] = min(min_e_[i], cost_function(e)[i]);
        }
    }
    
    NodeArray<NodeEntry> node_entries(graph);
    
    list<node> queue;
    
    {
        Label initial_label(Point(0.0, dimension),
                            source,
                            *this);
        
        node_entries[source].push_back(std::move(initial_label));
    }
    
    queue.push_back(source);
    
    while(!queue.empty()) {
        node current_node = queue.front();
        queue.pop_front();
        
        NodeEntry& current_node_entry = node_entries[current_node];
        
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
                                    *this);
                    
                    new_labels.push_back(std::move(new_label));
                }
                
                
                bool changed = check_domination(new_labels,
                                                neighbor_entry);
                
                if(changed) {
                    queue.push_back(neighbor);
                }
                
                
            }

            current_node_entry.proceed_labels_it();
        }
        
    }
    
    for(auto label : node_entries[target].labels()) {
        std::cout << label.cost << std::endl;
    }
    
    std::cout << "no labels: " << node_entries[target].labels().size() << std::endl;
}

}