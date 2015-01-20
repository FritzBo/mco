//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "est_ok10_module.h"

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

#include <mco/benchmarks/kurz_parser.h>
#include <mco/basic/point.h>
#include <mco/est/ok10/est_ok10_solver.h>

using mco::Point;
using mco::KurzParser;
using mco::EstOk10Solver;

void EstOk10Module::perform(int argc, char** argv) {
    try {
        CmdLine cmd("OK10 method to find the XND part of the Pareto-frontier of the efficient spanning tree problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-8, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();

        Graph graph;
        EdgeArray<Point *> edge_array(graph);

        KurzParser parser(file_name);
        parser.get_graph(graph, edge_array);

        assert(graph.numberOfEdges() > 0);

        unsigned objectives = edge_array[graph.chooseEdge()]->dimension();

        EstOk10Solver solver(epsilon);

        solver.Solve(graph, edge_array, objectives);

        solutions_.insert(solutions_.begin(),
                          solver.solutions().cbegin(),
                          solver.solutions().cend());

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EstOk10Module::solutions() {
    return solutions_;
}

string EstOk10Module::statistics() {
    string stats("");
    return stats;
}
