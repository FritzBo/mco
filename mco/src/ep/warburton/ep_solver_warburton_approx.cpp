/*
 * warburton.cpp
 *
 *  Created on: 02.04.2013
 *      Author: fritz
 */

#include <vector>
#include <limits>
#include <cmath>
#include <functional>
#include <list>
#include <iostream>

using std::vector;
using std::min;
using std::max;
using std::log;
using std::pow;
using std::ceil;
using std::floor;
using std::numeric_limits;
using std::function;
using std::list;
using std::pair;
using std::cout;
using std::endl;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/basic/ep_lagrange_bisect.h>
#include <mco/ep/martins/martins.h>
#include <mco/ep/warburton/ep_solver_warburton_approx.h>

namespace mco {

bool EpSolverWarburtonApprox::Solve(ogdf::Graph& graph,
                                    std::function<Point&(edge)> cost_function,
                                    unsigned dimension,
                                    const ogdf::node source,
                                    const ogdf::node target,
                                    const Point &epsilon,
                                    const Point& bound,
                                    bool test_only,
                                    bool directed,
                                    unsigned int processes,
                                    double theta) {
    
	Point ub(numeric_limits<double>::infinity(), dimension);
	Point lb(-numeric_limits<double>::infinity(), dimension);
    
	Point label_limits(numeric_limits<double>::infinity(), dimension);

	// Computing the bounds
    if(!computeBounds(graph,
                      cost_function,
                      epsilon,
                      bound,
                      source,
                      target,
                      theta,
                      dimension,
                      directed,
                      lb,
                      ub,
                      label_limits)) {
        
        return false;
    }

    const unsigned number_nodes = graph.numberOfNodes();
    unsigned skip_function = compute_skip_function(dimension, lb, ub);
    label_limits[skip_function] = bound[skip_function];
    
    cout << label_limits << endl;
    
    if(test_only) {
        unsigned number_subproblems = 1;
        
        cout << "bounds:" << endl;
        for(unsigned int k = 0; k < dimension; ++k) {
            cout << "k: " << k << ", lb: " << lb[k] << ", ub: " << ub[k] << ", limit: " << label_limits[k] << endl;
            
            if(k != skip_function) {
                number_subproblems *= static_cast<unsigned>(round(max(1.0, ub[k] - lb[k])));
            }
            
        }
        
        cout << "Number of subproblems to be solved: " << number_subproblems << endl;
        
        return true;
    } else {
#ifndef NDEBUG
	cout << "bounds:" << endl;
	for(unsigned int k = 0; k < dimension - 1; ++k)
		cout << "k: " << k << ", lb: " << lb[k] << ", ub: " << ub[k] << ", limit: " << label_limits[k] << endl;
#endif
    }
    
    unsigned k = skip_function == 0 ? 1 : 0;
    Point current_log_bound(lb);
    
    current_log_bound[skip_function] = numeric_limits<double>::quiet_NaN();
    
    list<pair<list<edge>, Point>> solutions;
    
    list<Point> nondominated_cells;
    
    Dijkstra<double, PairComparator<double, std::less<double>>> sssp_solver;
    NodeArray<double> distances(graph);
	NodeArray<edge> predecessor(graph);
    
    while(true) {
        
        bool dominated_cell = false;
        for(auto& cell : nondominated_cells) {
            bool dominates = true;
            for(unsigned i = 0; i < dimension; ++i) {
                if(i != skip_function) {
                    if(cell[i] >= current_log_bound[i]) {
                        dominates = false;
                        break;
                    }
                }
            }
            if(dominates) {
                dominated_cell = true;
                break;
            }
        }
        
        if(!dominated_cell) {
            
//#ifndef NDEBUG
            cout << "Computing for bound: " << current_log_bound << endl;
//#endif

        
            list<Point> new_bounds;
            Point T(dimension);
            for(unsigned i = 0; i < dimension; ++i) {
                if(i != skip_function) {
                    T[i] = epsilon[i] * pow(theta, current_log_bound[i]) / (number_nodes - 1);
                    
                }
                
                Point new_bound(0.0, dimension + 1);
                new_bound[i] = 1.0;
                new_bound[dimension] = - label_limits[i];
                new_bounds.push_back(std::move(new_bound));
            }
            
            EdgeArray<Point> scaled_costs(graph, Point(dimension));
            
            for(auto e : graph.edges) {
                Point scaled_cost(cost_function(e));
                for(unsigned i = 0; i < dimension; ++i) {
                    if(i != skip_function) {
                        scaled_cost[i] = floor(cost_function(e)[i] / T[i]);
                    }
                }

                scaled_costs[e] = std::move(scaled_cost);
            }
            
            auto scaled_cost_function = [&scaled_costs] (edge e) {
                return &scaled_costs[e];
            };
            
            vector<NodeArray<double>> heuristic_lower_bounds(dimension, graph);
            
            for(unsigned k = 0; k < dimension; ++k) {
                
                auto weight = [scaled_cost_function, k] (edge e) {
                    return scaled_cost_function(e)->operator[](k);
                };
                
                sssp_solver.singleSourceShortestPaths(graph,
                                                      weight,
                                                      target,
                                                      predecessor,
                                                      heuristic_lower_bounds[k],
                                                      directed ? DijkstraModes::Backward :
                                                      DijkstraModes::Undirected);
            }
            
            
            
            auto heuristic = [&heuristic_lower_bounds] (node n, unsigned objective) {
                return heuristic_lower_bounds[objective][n];
            };
            
            
            if(!lagrange_prune(graph,
                               scaled_cost_function,
                               source,
                               target,
                               directed,
                               label_limits,
                               skip_function,
                               dimension
               )) {
            
                EpSolverMartins solver;
                
                solver.Solve(graph,
                             scaled_cost_function,
                             dimension,
                             source,
                             target,
                             Point(numeric_limits<double>::infinity(), dimension),
                             std::list<std::pair<ogdf::NodeArray<Point*>,
                             ogdf::NodeArray<ogdf::edge>>>(),
                             heuristic,
                             new_bounds,
                             directed);
                
                unsigned new_solution_count = 0;
                bool nondominated_cell = false;
                
                for(auto solution_pair : solver.solutions()) {
                    auto solution = solution_pair.first;
                    
                    Point value(dimension);
                    for(auto e : solution) {
                        value += cost_function(e);
                    }

#ifndef NDEBUG
                    cout << "New solution: " << value << endl;
#endif
                
                    bool dominated = false;
                    auto sol_it = solutions.begin();
                    while(sol_it != solutions.end()) {
                        if(ComponentwisePointComparator(0, false)(sol_it->second, value)) {
                            dominated = true;
                            break;
                        }
                        
                        if(ComponentwisePointComparator(0, false)(value, sol_it->second)) {
                            sol_it = solutions.erase(sol_it);
                        }
                        
                        sol_it++;
                    }
                    if(!dominated) {
                        solutions.push_back(make_pair(solution, value));
                        nondominated_cell = true;
                    }
                }
                
                if(nondominated_cell) {
                    nondominated_cells.push_back(current_log_bound);
                }
            
#ifndef NDEBUG
                cout << "Number of new solutions acquired: " <<
                    new_solution_count << endl;
#endif
            
            }  else {
#ifndef NDEBUG
                cout << "Lagrange pruned." << endl;
#endif
            }
            
        }
        
        if(current_log_bound[k] < ub[k] - 1) {
            current_log_bound[k] += 1;
        } else {
            
            while(k < dimension && current_log_bound[k] >= ub[k] - 1) {
                k++;
                if(k == skip_function) {
                    k++;
                }
            }
            
            if(k == dimension) {
                break;
            }
            
            current_log_bound[k] += 1;
            
            for(unsigned i = 0; i < k; ++i) {
                if(i != skip_function) {
                    current_log_bound[i] = lb[i];
                }
            }
            
            k = skip_function == 0 ? 1 : 0;
        }
    }
    
    add_solutions(solutions.begin(), solutions.end());
    
    return true;

}
    
unsigned EpSolverWarburtonApprox::
compute_skip_function(unsigned dimension,
                      const Point& lb,
                      const Point& ub) {
    unsigned skip_function = 0;
    int no_trials = 0;
    for(unsigned i = 0; i < dimension; ++i) {
        if(max((double) no_trials, ub[i] - lb[i]) > no_trials) {
            no_trials = ub[i] - lb[i];
            skip_function = i;
        }
    }
    return skip_function;
}
    
bool EpSolverWarburtonApprox::
lagrange_prune(const Graph& graph,
               function<Point*(edge)> cost_function,
               const node source,
               const node target,
               bool directed,
               Point label_limits,
               unsigned skip_function,
               unsigned dimension) {
    
    EpLagrangeBisect bisect;
    
    Point* lambda = nullptr;
    Point* temp_label_limits = nullptr;
    
    function<Point&(edge)> lagrange_cost_function;
    
    std::list<Point> points;
    
    if(label_limits[skip_function] == numeric_limits<double>::infinity()) {
        
        lagrange_cost_function = [cost_function, skip_function, dimension, &points] (edge e) -> Point&{
            points.emplace_front(dimension - 1);
            Point& new_point = *points.begin();
            Point& target_point = *cost_function(e);
            unsigned j = 0;
            for(unsigned i = 0; i < dimension; ++i) {
                if(i != skip_function) {
                    new_point[j] = target_point[i];
                    j++;
                }
            }
            return *points.begin();
        };

        temp_label_limits = new Point(dimension - 1);
        unsigned j = 0;
        for(unsigned i = 0; i < dimension; ++i) {
            if(i != skip_function) {
                temp_label_limits->operator[](j) = label_limits[i];
                j++;
            }
        }
        
        lambda = new Point(dimension - 1);

        dimension--;
        
    } else {
        
        lagrange_cost_function = [cost_function] (edge e) -> Point& {
            return *cost_function(e);
        };
        
        lambda = new Point(dimension);
        temp_label_limits = new Point(label_limits);
        
    }
    
    double value = bisect.find_lagrange_multi(graph,
                                              lagrange_cost_function,
                                              dimension,
                                              source,
                                              target,
                                              directed,
                                              *temp_label_limits,
                                              *lambda);

    bool prune = (*temp_label_limits) * (*lambda) < value;
    
    delete lambda;
    delete temp_label_limits;
    
    return prune;
}
    
bool EpSolverWarburtonApprox::
computeBounds(const Graph& graph,
              std::function<Point&(edge)> cost_function,
              const Point& epsilon,
              const Point& bound,
              const node source,
              const node target,
              double theta,
              unsigned dimension,
              bool directed,
              Point& lb,
              Point& ub,
              Point& label_limits) {
    
    Point min_e(numeric_limits<double>::infinity(), dimension);
	Point max_e(-numeric_limits<double>::infinity(), dimension);
    
    double weight;
    
    const unsigned int number_nodes = graph.numberOfNodes();
    
    Dijkstra<double, PairComparator<double, std::less<double>>> sssp_solver;
    NodeArray<double> distances(graph);
	NodeArray<edge> predecessor(graph);

    for(unsigned int k = 0; k < dimension; ++k) {
        
#ifndef NDEBUG
        cout << "Objective function: " << k << endl;
#endif
        
        // If epsilon == 0 we have to divide by 0 and the approximation would
        // in fact be exact. So we disallow values of 0.
        if(epsilon[k] == 0) {
            return false;
        }
        
        // Let's find the lightest and heaviest edge
		for(auto e : graph.edges) {
			weight = cost_function(e)[k];
            if(weight != round(weight)) {
                return false;
            }
			min_e[k] = min(min_e[k], weight);
			max_e[k] = max(max_e[k], weight);
		}
        
#ifndef NDEBUG
        cout << "Heaviest edge: " << max_e[k] << endl;
        cout << "Lightest edge: " << min_e[k] << endl;
#endif
        
        // Initialize bounds for the labeling algorithm
		label_limits[k] = (number_nodes - 1) * theta / epsilon[k];
        
        // Compute an upper bound on the longest path
		ub[k] = min({ceil(log(max_e[k] * (number_nodes - 1)) / log(theta)),
            
            // Below this, all scaled edge weights in this objective
            // will be 0
            ceil(log(max_e[k] * (number_nodes - 1) / epsilon[k]) / log(theta)),
            
            // Using the bounds, the user gives us, but at least 1.0
            max(1.0, ceil(log(bound[k])/log(theta)))});
        
#ifndef NDEBUG
        cout << "Log-upperbound on the longest path: " << ub[k] << endl;
#endif
        
        // Since L_k > (n - 1) / \varepsilon (p. 75)
        lb[k] = max({1.0,
            floor(log((number_nodes - 1)/epsilon[k])/log(theta))});
        
        // Try to rise the lower bound to reduce the number of subproblems
        // to be solved.
        unsigned i = lb[k];
        while(lb[k] == i && ub[k] > lb[k]) {
            
#ifndef NDEBUG
            cout << "Lower bound to probe: " << i << endl;
#endif
            
            double Tk = epsilon[k] * pow(theta, i) / (number_nodes - 1);
            
            // Compute a scaled weight function for i
            auto scaled_weight_function = [cost_function, k, &Tk] (edge e) {
                return floor(cost_function(e)[k] / Tk);
            };
            
            // Compute shortest path for objective k and scaling of i
            sssp_solver.singleSourceShortestPaths(graph,
                                                  scaled_weight_function,
                                                  target,
                                                  predecessor,
                                                  distances,
                                                  directed ? DijkstraModes::Backward : DijkstraModes::Undirected);
            
#ifndef NDEBUG
            cout << "Shortest path distance: " << distances[source] << endl;
#endif
            
            // if the shortest path with respect to objective k is longer
            // than the absolute limit on objective k: rise lower bound by one.
            // Else: We found a lower bound which results in a feasible
            // instance.
            if(distances[source] > label_limits[k]) {
                lb[k]++;
            }
            
            i++;
        }
        
#ifndef NDEBUG
        cout << "===============================" << endl;
#endif
        
	}
    return true;
}

} // namespace mco
