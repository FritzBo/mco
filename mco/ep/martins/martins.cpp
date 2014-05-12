/*
 * martins.cpp
 *
 *  Created on: 23.03.2013
 *      Author: fritz
 */

#include <mco/ep/martins/martins.h>

#include <queue>
#include <vector>
#include <set>
#include <list>

using std::priority_queue;
using std::vector;
using std::set;
using std::list;
using std::function;
using std::pair;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::AdjElement;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;

#include <mco/basic/point.h>
#include <mco/ep/basic/ep_instance.h>
#include <mco/ep/martins/label.h>

namespace mco {
    
void EpSolverMartins::
Solve(Graph& graph,
      function<const Point*(edge)> weights,
      unsigned dimension,
      node source,
      node target,
      const Point& bound,
      function<double(ogdf::node, unsigned)> heuristic,
      list<pair<NodeArray<Point*>,
                NodeArray<edge>>> initial_labels,
      list<Point> first_phase_bounds,
      bool directed) {
    
    using LabelPriorityQueue = priority_queue<Label *, vector<Label *>, LexLabelComp>;
    
    unsigned bound_deletion = 0;
    unsigned heuristic_deletion = 0;
    unsigned first_phase_deletion = 0;
    
	LabelPriorityQueue lex_min_label((LexLabelComp()));
	NodeArray<list<Label *>> labels(graph);
    
    ComponentwisePointComparator comp_leq(epsilon_, false);

	Label *null_label = new Label(Point::Null(dimension), source, nullptr);
    null_label->in_queue = true;
	labels[source].push_back(null_label);
    
	lex_min_label.push(null_label);
    
    if(!initial_labels.empty()) {
        construct_labels(labels, initial_labels);
        
        for(auto n : graph.nodes) {
            if(n != target && n != source) {
                for(auto label : labels[n]) {
                    lex_min_label.push(label);
                }
            }
        }
    }

	while(!lex_min_label.empty()) {
		Label *label = lex_min_label.top();
		lex_min_label.pop();
        assert(label->in_queue);
        label->in_queue = false;

		if(label->mark_dominated) {
			delete label;
			continue;
		}

		const Point *label_cost = label->point;
		node n = label->n;
        
        if(n == target) {
//            cout << *label_cost << endl;
            continue;
        }

//		cout << endl << n << ", " << *label->point << ": ";

		AdjElement *adj;
		forall_adj(adj, n) {
			edge e = adj->theEdge();

			if(e->isSelfLoop())
				continue;

			node v = e->target();

            if(directed) {
                if(v == n)
                    continue;
            } else {
                if(v == n) {
                    v = e->source();
                }
            }

			if(v == source)
				continue;

//			cout << v << ", ";

			const Point *edge_cost = weights(e);
			const Point *new_cost = new Point(*edge_cost + *label_cost);			// Owner
            
            for(unsigned i = 0; i < dimension; ++i) {
                double heuristic_cost = new_cost->operator[](i) + heuristic(v, i);
                if(heuristic_cost > bound[i]) {
                    delete new_cost;
                    new_cost = nullptr;
                    ++bound_deletion;
                    break;
                }
                
            }
            
            if(new_cost == nullptr) {
                continue;
            }
            
            for(auto label : labels[target]) {
                const Point& cost = *label->point;
                
                bool dominated = true;
                for(unsigned i = 0; i < dimension; ++i) {
                    if(new_cost->operator[](i) + heuristic(v, i) < cost[i]) {
                        dominated = false;
                        break;
                    }
                }
                
                if(dominated) {
                    delete new_cost;
                    new_cost = nullptr;
                    ++heuristic_deletion;
                    break;
                }
            }
            
            if(new_cost == nullptr) {
                continue;
            }
            
            for(auto cost : first_phase_bounds) {
                bool dominated = true;
                for(unsigned i = 0; i < dimension; ++i) {
                    if(new_cost->operator[](i) + heuristic(v, i) <= cost[i]) {
                        dominated = false;
                        break;
                    }
                }
                
                if(dominated) {
                    delete new_cost;
                    new_cost = nullptr;
                    ++first_phase_deletion;
                    break;
                }
            }
            
            if(new_cost == nullptr) {
                continue;
            }

			bool dominated = false;

			auto iter = labels[v].begin();
			while(iter != labels[v].end()) {

				Label * target_label = *iter;

				if(comp_leq(target_label->point, new_cost)) {
					dominated = true;
					break;
				}

                if(target_label->in_queue) {
                    if(comp_leq(new_cost, target_label->point)) {
                        target_label->mark_dominated = true;
                        iter = labels[v].erase(iter);
                    } else {
                        ++iter;
                    }
                } else {
                    ++iter;
                }
			}

			if(dominated) {
				delete new_cost;
				continue;
			}

			Label * new_label = new Label(new_cost, v, label);
			labels[v].push_back(new_label);

			lex_min_label.push(new_label);
            new_label->in_queue = true;
		}

	}

	list<const Point *> target_labels;
	for(auto &label : labels[target])
		target_labels.push_back(new Point(*label->point));

	reset_solutions();
	add_solutions(target_labels.begin(), target_labels.end());

	for(Label *label : labels[target]) {
		const Label *current_label = label;
//		cout << "(";
		while(current_label->n != source) {
			cout << current_label->n << ", ";
			current_label = current_label->pred;
		}
		cout << current_label->n;
//        cout << ")";
        cout << endl;
        cout << *label->point << endl;
	}
    
    cout << "Length bound deletions: " << bound_deletion << endl;
    cout << "Heuristic bound deletions: " << heuristic_deletion << endl;
    cout << "First phase bound deletions: " << first_phase_deletion << endl;
    
	node n;
	forall_nodes(n, graph) {
		for(auto &label : labels[n])
			delete label;
	}
}
    
void EpSolverMartins::
construct_labels(NodeArray<list<Label*>> & labels,
                 list<pair<NodeArray<Point*>,
                           NodeArray<edge>>>& initial_labels) {
    
    LexPointComparator comp;
    EqualityPointComparator eq;
    
    for(auto solution : initial_labels) {
        auto distances = solution.first;
        auto predecessor = solution.second;
        NodeArray<bool> in_queue(*labels.graphOf(), true);
        
        auto order = [distances, comp] (node v, node w) {
            return comp(distances[v], distances[w]);
        };
        
        priority_queue<node, vector<node>, decltype(order)> nodes(order);
        
        for(auto n : labels.graphOf()->nodes) {
            nodes.push(n);
        }
        
        // For each node, starting with one with a highest label
        while(!nodes.empty()) {
            node n = nodes.top();
            nodes.pop();
            
            // If the node was already marked to be removed from the queue,
            // draw next node,
            // else, mark it as removed and proceed
            if(!in_queue[n]) {
                continue;
            }
            
            bool labeling_finished = false;
            
            list<node> path;
            Label * pred = nullptr;
            
            // While we are still scanning for the first label for the path and
            // we have not run into the end of the path and there are still
            // nodes on the path which are in the queue
            while(!labeling_finished || in_queue[n]) {
                
//                cout << n << endl;
                
                // If we are still looking for the root of the path
                if(!labeling_finished) {
                    // Is there a label with the same point?
                    for(auto label : labels[n]) {
                        if(eq(label->point, distances[n])) {
                            labeling_finished = true;
                            pred = label;
                        }
                    }
                    
                    // No label with the same point? Add the node
                    // to the path
                    if(!labeling_finished) {
                        path.push_back(n);
                    }
                }
                
                // Remove node from queue
                in_queue[n] = false;
                
                if(predecessor[n] == nullptr) {
                    break;
                }
                
                // Proceed with predecessor on current path
                n = n == predecessor[n]->source() ? predecessor[n]->target() :
                            predecessor[n]->source();
            }
            
            // For each node on the path
            while(!path.empty()) {
                // beginning with the latest added
                auto n = path.back();
                path.pop_back();
                
                // Create a new label
                auto label = new Label(distances[n], n, pred);
                // add it to the node
                labels[n].push_back(label);
                
                // And set the predecessor accordingly
                pred = label;
            }
            
        }
        
    }
}


}
