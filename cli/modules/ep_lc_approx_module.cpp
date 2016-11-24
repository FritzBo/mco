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
using std::vector;
using std::pair;
using std::function;
using std::numeric_limits;
using std::cout;
using std::endl;

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
using TCLAP::MultiArg;

#include <mco/ep/basic/instance_scalarizer.h>
#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/preprocessing/constrained_reach.h>
#include <mco/ep/lc_approx/lc_approx.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>
#include <mco/cli/parse_util.h>

using mco::LCApprox;
using mco::TemporaryGraphParser;
using mco::Point;
using mco::Dijkstra;
using mco::DijkstraModes;
using mco::ConstrainedReachabilityPreprocessing;
using mco::InstanceScalarizer;

void EpLCApproxModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Label Correcting Approximation for the Efficient Path Problem.", ' ', "0.1");
        
        MultiArg<string> epsilon_argument("E", "epsilon", "The approximation factor to use.", false, "epsilon");
        
        ValueArg<double> epsilon_all_argument("e", "epsilon-all", "The approximation factor to use for all objective functions.", false, 0.0, "epsilon");

        ValueArg<unsigned> exact_argument("X", "exact", "When using the epsilon-all or e parameter, this objective function will not be approximated.", false, 0, "objective function");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);

        cmd.add(epsilon_argument);
        cmd.add(epsilon_all_argument);
        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(exact_argument);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        bool directed = is_directed_arg.getValue();
        double epsilon_all = epsilon_all_argument.getValue();
        unsigned exact_objective = exact_argument.getValue();
        
        Graph graph;
        EdgeArray<Point> raw_costs(graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, graph, raw_costs, dimension, source, target);

        EdgeArray<Point> costs(graph, Point(dimension));
        Point factor(100.0, dimension);
        InstanceScalarizer::scaleround_instance(graph,
                                                raw_costs,
                                                dimension,
                                                factor,
                                                costs);


        vector<NodeArray<double>> distances(dimension, graph);

        auto cost_function = [&costs] (edge e) -> Point& {
            return costs[e];
        };
        
        Point epsilon(epsilon_all, dimension);
        if(exact_argument.isSet()) {
            epsilon[exact_objective - 1] = 0.0;
        }
        
        parse_epsilon(epsilon_argument,
                      dimension,
                      epsilon);
        
        LCApprox solver;

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

void EpLCApproxModule::parse_epsilon(const MultiArg<string>& epsilon_argument,
                                     unsigned dimension,
                                     Point& epsilon) {
    
    auto epsilon_it = epsilon_argument.begin();
    while(epsilon_it != epsilon_argument.end()) {
        vector<string> tokens;
        mco::tokenize(*epsilon_it, tokens, ":");
        
        if(tokens.size() != 2) {
            cout << "Error" << endl;
            return;
        }
        
        unsigned objective_function = stoul(tokens[0]);
        double epsilon_value = stod(tokens[1]);
        
        if(objective_function > dimension) {
            cout << "Error" << endl;
            return;
        }
        
        if(objective_function == 0) {
            cout << "Error" << endl;
            return;
        }
        
        epsilon[objective_function - 1] = epsilon_value;
        
        epsilon_it++;
    }
    
}

const list<pair<const list<edge>, const Point>>& EpLCApproxModule::solutions() {
    return solutions_;
}

string EpLCApproxModule::statistics() {
    string stats("");
    return stats;
}
