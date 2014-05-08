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

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::list;

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
    
    auto cost_function = [costs] (edge e) {
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
        
        for(auto p : solver.solutions()) {
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
            auto length = [costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [distances] (node n, unsigned objective) {
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
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }

    } else if(algorithm.compare("fpre-martins") == 0) {
        Dijkstra<double> sssp_solver;
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<edge> predecessor(graph);
        
        cout << "calculating heuristic..." << endl;
        
        for(unsigned i = 0; i < dimension; ++i) {
            auto length = [costs, i] (edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  length,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  DijkstraModes::Undirected);
        }
        
        auto heuristic = [distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        Point bound(dimension);
        for(unsigned i = 0; i < dimension - 1; ++i) {
            bound[i] = numeric_limits<double>::infinity();
        }
        bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
        
        cout << "Running first phase..." << endl;
        
        EPDualBensonSolver<> weighted_solver;
        
        weighted_solver.Solve(graph, cost_function, source, target);
        
        list<Point> first_phase_bounds;
        
        for(auto point : weighted_solver.solutions()) {
            first_phase_bounds.push_back(*point);
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
        
        cout << solver.solutions().size() << endl;
        
        for(auto p : solver.solutions()) {
            cout << *p << endl;
        }
        
    } else {
        cout << "Unknown algorithm: " << algorithm << endl;
    }
    
    
}