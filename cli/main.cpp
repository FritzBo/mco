//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 09.04.14.
//
//

#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::node;
using ogdf::edge;

#include <mco/basic/point.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/ep/brum_shier/ep_solver_bs.h>
#include <mco/ep/brum_shier/ep_weighted_bs.h>
#include <mco/ep/brum_shier/ep_weighted_bs.h>
#include <mco/ep/martins/martins.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/generic/benson_dual/ove_node_lists.h>

using mco::Point;
using mco::TemporaryGraphParser;
using mco::EpSolverBS;
using mco::EpWeightedBS;
using mco::EPDualBensonSolver;
using mco::EpSolverMartins;

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
        cout << "Not yet implemented" << endl;
    } else {
        cout << "Unknown algorithm: " << algorithm << endl;
    }
    
    
}