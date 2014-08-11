//
//  ep_lc_approx_module.cpp
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#include "ep_lc_approx_module.h"

using std::map;
using std::string;
using std::list;
using std::pair;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;

#include <mco/ep/lc_approx/lc_approx.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>

using mco::LCApprox;
using mco::TemporaryGraphParser;
using mco::Point;

void EpLCApproxModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Label Correcting Approximation for the Efficient Path Problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "The approximation factor to use.", true, 1, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);
        
        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        
        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();
        bool directed = is_directed_arg.getValue();
        
        Graph graph;
        EdgeArray<Point> costs(graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, graph, costs, dimension, source, target);
        
        LCApprox solver;
        
        auto cost_function = [&costs] (edge e) -> Point& {
            return costs[e];
        };
        
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     directed,
                     epsilon);
        
        solutions_.insert(solutions_.begin(),
                          solver.solutions().cbegin(),
                          solver.solutions().cend());
        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EpLCApproxModule::solutions() {
    return solutions_;
}

string EpLCApproxModule::statistics() {
    string stats("");
    return stats;
}
