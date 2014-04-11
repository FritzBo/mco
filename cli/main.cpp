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

using mco::Point;
using mco::TemporaryGraphParser;
using mco::EpSolverBS;

int main(int argc, char** argv) {
//    if(argc != 2) {
//        cout << "Usage: " << argv[0] << " <file>" << endl;
//    }
    
//    string filename(argv[1]);
    string filename("../../instances/ep/graph_1000_1000");
//    string filename("/Users/fritz/Desktop/grid70_100_7");
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
    
    EpSolverBS solver;
    
    auto cost_function = [costs] (edge e) {
        return &costs[e];
    };
    
    solver.Solve(graph, cost_function, dimension, source, target, false);
    
    cout << solver.solutions().size() << endl;
    
    for(auto p : solver.solutions()) {
        cout << *p << endl;
    }
    
    
}