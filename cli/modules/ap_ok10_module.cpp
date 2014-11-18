//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "ap_ok10_module.h"

#include <set>

using std::set;
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

#include <mco/benchmarks/mcap_parser.h>
#include <mco/ap/ok10/ap_ok10_solver.h>
#include <mco/basic/point.h>

using mco::Point;
using mco::MCAPParser;
using mco::AssignmentInstance;
using mco::ApOk10Solver;

void ApOk10Module::perform(int argc, char** argv) {
    try {
        CmdLine cmd("OK10 method to find the XND part of the Pareto-frontier of the efficient assignment problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-8, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();

        Graph graph;
        EdgeArray<Point *> edge_array(graph);
        set<node> agents;

        MCAPParser parser(file_name);
        AssignmentInstance instance = parser.get_instance(graph,
                                                          edge_array,
                                                          agents);

        ApOk10Solver solver(epsilon);

        solver.Solve(instance);

        solutions_.insert(solutions_.begin(),
                          solver.solutions().cbegin(),
                          solver.solutions().cend());

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& ApOk10Module::solutions() {
    return solutions_;
}

string ApOk10Module::statistics() {
    string stats("");
    return stats;
}
