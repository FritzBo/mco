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

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <mco/ep/basic/dijkstra.h>
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
    
	const unsigned int number_nodes = graph.numberOfNodes();

//	vector<EdgeArray<double>> weight_functions(dimension - 1, EdgeArray<double>(graph, 0));
	NodeArray<double> distances(graph);
	NodeArray<edge> predecessor(graph);
    
	Point min_e(numeric_limits<double>::infinity(), dimension);
	Point max_e(-numeric_limits<double>::infinity(), dimension);
    
	Point ub(numeric_limits<double>::infinity(), dimension);
	Point lb(-numeric_limits<double>::infinity(), dimension);
    
	Point label_limits(numeric_limits<double>::infinity(), dimension);
    
	double weight;

    Dijkstra<double> sssp_solver;

	// Computing the bounds
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
    
    unsigned skip_function = 0;
    int no_trials = 0;
    for(unsigned i = 0; i < dimension; ++i) {
        if(max((double) no_trials, ub[i] - lb[i]) > no_trials) {
            no_trials = ub[i] - lb[i];
            skip_function = i;
        }
    }
    
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
                    
                    Point new_bound(0.0, dimension + 1);
                    new_bound[i] = 1.0;
                    new_bound[dimension] = - label_limits[i];
                    new_bounds.push_back(std::move(new_bound));
                }
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
            
            k = 0;
        }
    }
    
    add_solutions(solutions.begin(), solutions.end());
    
    return true;

}

} // namespace mco
