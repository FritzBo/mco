//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 09.04.14.
//
//

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <chrono>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::list;
using std::pair;
using std::make_pair;
using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <mco/basic/point.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/brum_shier/ep_solver_bs.h>
#include <mco/ep/brum_shier/ep_weighted_bs.h>
#include <mco/ep/brum_shier/ep_weighted_bs.h>
#include <mco/ep/martins/martins.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/ep/martins/weighted_martins.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/generic/benson_dual/ove_node_lists.h>

using mco::Point;
using mco::TemporaryGraphParser;
using mco::EpSolverBS;
using mco::EpWeightedBS;
using mco::EPDualBensonSolver;
using mco::EpSolverMartins;
using mco::EpWeightedMartins;
using mco::Dijkstra;
using mco::DijkstraModes;
using mco::LexPointComparator;
using mco::ParetoDominationPointComparator;
using mco::ComponentwisePointComparator;
using mco::EqualityPointComparator;

int main(int argc, char** argv) {
    if(argc != 3) {
        cout << "Usage: " << argv[0] << "<algorithm> <file>" << endl;
    }
    
    string algorithm(argv[1]);
    string filename(argv[2]);

    Graph graph;
    EdgeArray<Point> costs(graph);
    node source;
    node target;
    unsigned dimension;
    
    mco::TemporaryGraphParser parser;
    
    parser.getGraph(filename,
                    graph,
                    costs,
                    dimension,
                    source,
                    target);
    
    cout << "Read a graph with " << graph.numberOfNodes() <<
    " nodes and " << graph.numberOfEdges() << " edges." << endl;
    cout << "Solving..." << endl;
    
    auto cost_function = [&costs] (edge e) {
        return &costs[e];
    };
    
    if (algorithm.compare("dual-benson") == 0) {
        
        EPDualBensonSolver<> solver;
        
        solver.Solve(graph, cost_function, source, target);
        
        cout << solver.solutions().size() << endl;
    
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
        
    } else if(algorithm.compare("label-correcting") == 0) {
    
        EpSolverBS solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        vector<const Point*> solutions(solver.solutions().cbegin(),
                                       solver.solutions().cend());
        
        std::sort(solutions.begin(), solutions.end(), mco::LexPointComparator());
        
//        mco::ComponentwisePointComparator comp;
//        mco::EqualityPointComparator eq;
//        bool error = false;
//        
//        auto iter1 = solutions.begin();
//        while(iter1 != solutions.end()) {
//            auto iter2 = iter1 + 1;
//            while(iter2 != solutions.end()) {
//                if(comp(*iter1, *iter2) || comp(*iter2, *iter1) || eq(*iter1, *iter2)) {
//                    cout << *(*iter1) << endl;
//                    cout << *(*iter2) << endl;
//                    error = true;
//                }
//                ++iter2;
//            }
//            ++iter1;
//        }
//        
//        cout << error << endl;
        
        for(auto p : solutions) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("w-label-correcting") == 0) {
        EpWeightedBS solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("martins") == 0) {
        EpSolverMartins solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("w-martins") == 0) {
        EpWeightedMartins solver;
        
        solver.Solve(graph, cost_function, dimension, source, target, false);
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
    } else if(algorithm.compare("pre-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "calculating heuristic..." << endl;
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [&costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [&distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        Point bound(dimension);
        for(unsigned i = 0; i < dimension - 1; ++i) {
            bound[i] = numeric_limits<double>::infinity();
        }
        bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        EpSolverMartins solver;
        
        cout << "Running Martins algorithm..." << endl;
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     bound,
                     heuristic,
                     list<Point>(),
                     false);
        
        cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
        
    } else if(algorithm.compare("fpre-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "calculating heuristic..." << endl;
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [&costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [&distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        Point bound(dimension);
        for(unsigned i = 0; i < dimension - 1; ++i) {
            bound[i] = numeric_limits<double>::infinity();
        }
        bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        cout << "Running first phase..." << endl;
        
        list<Point> first_phase_bounds;
        {
        
            EPDualBensonSolver<> weighted_solver;
            
            weighted_solver.Solve(graph, cost_function, source, target);
            
            
            for(auto point : weighted_solver.solutions()) {
                first_phase_bounds.push_back(*point);
            }
            
        }
        
        EpSolverMartins solver;
        
        cout << "Running Martins algorithm..." << endl;
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     bound,
                     heuristic,
                     first_phase_bounds,
                     false);
        
        cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
        
    } else if(algorithm.compare("flabel-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "Calculating heuristic... ";
        steady_clock::time_point start = steady_clock::now();
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [&costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [&distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        steady_clock::time_point heuristic_end = steady_clock::now();
        duration<double> heuristic_computation_span
            = duration_cast<duration<double>>(heuristic_end - start);
        
        cout << "Done. (" << heuristic_computation_span.count() << "s)" << endl;
        
        Point bound(dimension);
        for(unsigned i = 0; i < dimension; ++i) {
            bound[i] = numeric_limits<double>::infinity();
        }
//        bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        cout << "Running first phase... ";
        
        list<pair<NodeArray<Point *>, NodeArray<edge>>> solutions;
        
        auto callback = [&solutions, &graph, dimension] (NodeArray<Point *>& distances,
                                     NodeArray<edge>& predecessors) {
            
            NodeArray<Point *> new_distances(graph);
            
            for(auto n : distances.graphOf()->nodes) {
                Point* p = new Point(dimension);
                std::copy(distances[n]->cbegin() + 1, distances[n]->cend(),
                          p->begin());
                new_distances[n] = p;
            }
            
            solutions.push_back(make_pair(new_distances, predecessors));
        };
        
        {
            
            EPDualBensonSolver<> weighted_solver;
            
            weighted_solver.Solve(graph, cost_function, source, target, callback);
            
        }
        
        steady_clock::time_point fp_end = steady_clock::now();
        duration<double> fp_computation_span
            = duration_cast<duration<double>>(fp_end - heuristic_end);
        
        cout << "Done. (" << fp_computation_span.count() << "s)" << endl;

        vector<list<node>> paths;
        vector<Point> values;
        
        auto path_callback = [&paths] (list<node> path) {
            paths.push_back(path);
        };
        
        auto value_callback = [&values] (Point p) {
            values.push_back(p);
        };
        
        EpSolverMartins solver;
        
        solver.set_path_callback(path_callback);
        solver.set_value_callback(value_callback);
        
        cout << "Running Martins algorithm... ";
        std::flush(cout);
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     bound,
                     solutions,
                     heuristic,
                     false);
        
        steady_clock::time_point martins_end = steady_clock::now();
        duration<double> martins_computation_span
            = duration_cast<duration<double>>(martins_end - fp_end);
        
        cout << "Done. (" << martins_computation_span.count() << "s)" << endl;

        cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
        
    } else {
        cout << "Unknown algorithm: " << algorithm << endl;
    }
    
    
}