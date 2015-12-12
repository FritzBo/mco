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

using ogdf::node;
using ogdf::edge;
using ogdf::Graph;
using ogdf::NodeArray;

namespace mco {

bool LCApprox::
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

        dominated = false;

        unsigned check_label_it = 0;
        while(check_label_it < neighbor_labels_end) {

            auto check_label = neighbor_labels[check_label_it];

            if(check_label->deleted) {

                neighbor_entry.erase(neighbor_labels,
                                     check_label_it);

                --neighbor_labels_end;

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

                --neighbor_labels_end;
                
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

                    --neighbor_new_labels_end;

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
            pred->succ_label.push_back(new_label);

            changed = true;
        }

        ++new_label_it;
    }
    
    return changed;
}
    
bool LCApprox::check_heuristic_prunable(const Label& label,
                                        const Point bound,
                                        const list<Point>& bounds) {
    
    Point heuristic_cost(dimension_);
    for(unsigned i = 0; i < dimension_; ++i) {
        heuristic_cost[i] = label.cost[i] + heuristic_(label.n, i);
    }
    
    compute_pos(heuristic_cost, heuristic_cost);

    ComponentwisePointComparator leq(0, false);

    if(!leq(heuristic_cost, bound)) {
        return true;
    }

    for(auto& bound : bounds) {
        if(leq(heuristic_cost, bound)) {
            return false;
        }
    }
    
    return bounds.size() > 0;
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
            min_e_[i] = std::min(min_e_[i], cost_function(e)[i]);
        }
    }

    list<Point> scaled_disj_bounds;
    Point scaled_bound(dimension_);

    if(use_bounds_) {
        compute_pos(bound_, scaled_bound);

        for(auto& bound : disj_bounds_) {
            Point scaled_disj_bound(dimension_);
            compute_pos(bound, scaled_disj_bound);
#ifndef NDEBUG
            std::cout << scaled_disj_bound << std::endl;
#endif
            scaled_disj_bounds.push_back(std::move(scaled_disj_bound));
        }
    }
    
    NodeArray<NodeEntry> node_entries(graph);
    
    list<node> queue;
    
    {
        auto initial_label = new Label(Point(0.0, dimension),
                                       source,
                                       nullptr,
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
            
            for(auto adj : current_node->adjEntries) {
                edge current_edge = adj->theEdge();

                if(current_edge->isSelfLoop()) {
                    continue;
                }

                if(directed && current_edge->target() == current_node) {
                    continue;
                }
                
                auto neighbor = current_edge->opposite(current_node);
                
                auto& neighbor_entry = node_entries[neighbor];

                vector<Label*> new_labels(current_node_entry.new_labels_end(), nullptr);

                unsigned size = 0;

                unsigned current_label_it = 0;
                while(current_label_it < current_node_entry.new_labels_end()) {

                    auto label = current_new_labels[current_label_it];

                    assert(label != nullptr);
                    assert(label->n == current_node);

                    if(!label->deleted) {

                        auto new_label = new Label(label->cost + cost_function(current_edge),
                                                   neighbor,
                                                   current_edge,
                                                   label,
                                                   *this);
                        
                        if(!use_bounds_ || !use_heuristic_ ||
                           !check_heuristic_prunable(*new_label,
                                                     scaled_bound,
                                                     scaled_disj_bounds)) {
                            
                            new_labels[size] = new_label;
                            ++size;
                        }
                    }

                    ++current_label_it;
                }
                
                if(size != 0) {
                
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

    for(unsigned i = 0; i < node_entries[target].new_labels_end(); ++i) {
        auto label = node_entries[target].new_labels()[i];

        list<edge> path;
        const Label* curr = label;
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