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
                                    list<Point> linear_bounds,
                                    unsigned int processes,
                                    double theta,
                                    bool directed) {
    
	const unsigned int number_nodes = graph.numberOfNodes();

//	vector<EdgeArray<double>> weight_functions(dimension - 1, EdgeArray<double>(graph, 0));
	NodeArray<double> distances(graph);
	NodeArray<edge> predecessor(graph);
    
	Point min_e(numeric_limits<double>::infinity(), dimension - 1);
	Point max_e(-numeric_limits<double>::infinity(), dimension - 1);
    
	Point ub(numeric_limits<double>::infinity(), dimension - 1);
	Point lb(-numeric_limits<double>::infinity(), dimension - 1);
    
	Point label_limits(numeric_limits<double>::infinity(), dimension);
    
	double weight, d;
	int max_i, min_i;
    
    Dijkstra<double> sssp_solver;

	// Computing the bounds
	for(unsigned int k = 0; k < dimension - 1; ++k) {
        
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
			min_e[k] = min(min_e[k], weight);
			max_e[k] = max(max_e[k], weight);
		}
      
#ifndef NDEBUG
        cout << "Heaviest edge: " << max_e[k] << endl;
        cout << "Lightest edge: " << min_e[k] << endl;
#endif

        // Pre-initialize limits for labeling algorithm
		label_limits[k] = static_cast<int>(ceil((number_nodes - 1) * theta / epsilon[k]));

        // Compute an upper bound on the longest path
		max_i = static_cast<int>(ceil(
                                      log(
                                          max( 1.0,
                                              min(
                                                  max_e[k] * (number_nodes - 1),
                                                  max_e[k] * (number_nodes - 1) / epsilon[k])
                                              ) / log(theta)
                                          )
                                      ) + 1);
      
#ifndef NDEBUG
        cout << "Log-upperbound on the longest path: " << max_i << endl;
#endif
        
        // Compute a upper bound on the shortest path for generic epsilon
        min_i = static_cast<int>(log(
                                     min(
                                         static_cast<double>(max_i),
                                         max(
                                             {
                                                 1.0,
                                                 min_e[k] * (number_nodes - 1) / epsilon[k],
                                                 log(theta)
                                             })
                                         )
                                     )
                                 );
      
#ifndef NDEBUG
        cout << "Log-lowerbound on the shortest path: " << min_i << endl;
#endif

		for(int i = min_i; i < max_i; ++i) {
            
			d = static_cast<double>(pow(2, max_i - i));
      
#ifndef NDEBUG
            cout << "Probing bound: 2**" << (max_i - i) << endl;
#endif
            
            auto weight_function = [&cost_function, d, k, &epsilon, number_nodes] (edge e) {
                return static_cast<int>(floor(
                                              cost_function(e)[k] * (number_nodes - 1) / (epsilon[k] * d)
                                              ));
            };
            
			// TODO: Error Message / Exception
            sssp_solver.singleSourceShortestPaths(graph,
                                                  weight_function,
                                                  source,
                                                  predecessor,
                                                  distances,
                                                  directed ? DijkstraModes::Forward : DijkstraModes::Undirected);
            
            if(distances[target] == numeric_limits<double>::max()) {
                return false;
            }
      
#ifndef NDEBUG
            cout << "Shortest path: " << distances[target] << endl;
#endif

			if(distances[target] > (number_nodes - 1) * theta / epsilon[k] || max_i - i == 1) {
                
				lb[k] = max_i - i + 1;
				ub[k] = static_cast<int>(ceil(
                                              max(lb[k] - 1,
                                                  log(
                                                      min(
                                                          max_e[k] * (number_nodes - 1),
                                                          max_e[k] * (number_nodes - 1) / epsilon[k]) / log(theta)
                                                      )
                                                  ) + 1
                                              )
                                         );
				break;
                
			}
            
		}
        
#ifndef NDEBUG
        cout << "===============================" << endl;
#endif

	}
    
#ifndef NDEBUG
	cout << "bounds:" << endl;
	for(unsigned int k = 0; k < dimension - 1; ++k)
		cout << "k: " << k << ", lb: " << lb[k] << ", ub: " << ub[k] << ", limit: " << label_limits[k] << endl;
#endif
    
    unsigned k = 0;
    Point current_log_bound(lb);
    
    while(true) {
        
#ifndef NDEBUG
        cout << "Computing for bound: " << current_log_bound << endl;
#endif
        
        list<Point> new_bounds(linear_bounds);
        Point T(epsilon);
        for(unsigned i = 0; i < dimension - 1; ++i) {
            T[i] *= pow(theta, current_log_bound[i]) / (number_nodes - 1);
            
            Point new_bound(0.0, dimension + 1);
            new_bound[i] = 1.0;
            new_bound[dimension] = - label_limits[i];
            new_bounds.push_back(std::move(new_bound));
        }
        
        EdgeArray<Point> scaled_costs(graph, Point(dimension));
        
        for(auto e : graph.edges) {
            Point scaled_cost(cost_function(e));
            for(unsigned i = 0; i < dimension - 1; ++i) {
                scaled_cost[i] = floor(cost_function(e)[i] / T[i]);
            }
            scaled_cost[dimension - 1] = cost_function(e)[dimension - 1];
            scaled_costs[e] = std::move(scaled_cost);
        }
        
        auto scaled_cost_function = [&scaled_costs] (edge e) {
            return &scaled_costs[e];
        };
        
        auto heuristic = [] (node n, unsigned objective) {
            return 0.0;
        };
        
        EpSolverMartins solver;
        
        solver.Solve(graph,
                     scaled_cost_function,
                     dimension,
                     source,
                     target,
                     Point(numeric_limits<double>::infinity()),
                     std::list<std::pair<ogdf::NodeArray<Point*>,
                     ogdf::NodeArray<ogdf::edge>>>(),
                     heuristic,
                     new_bounds,
                     directed);
        
        for(auto solution_pair : solver.solutions()) {
            auto solution = solution_pair.first;
            auto value = solution_pair.second;
            
            cout << value << endl;
        }
        
        add_solutions(solver.solutions().begin(),
                      solver.solutions().end());
        
        if(current_log_bound[k] < ub[k] - 1) {
            current_log_bound[k] += 1;
        } else {
            
            while(k < dimension - 1 && current_log_bound[k] >= ub[k] - 1) {
                k++;
            }
            
            if(k == dimension - 1) {
                break;
            }
            
            current_log_bound[k] += 1;
            
            for(unsigned i = 0; i < k; ++i) {
                current_log_bound[i] = lb[i];
            }
            
            k = 0;
        }
    }
    
    return true;

}

} // namespace mco
