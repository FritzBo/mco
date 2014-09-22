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
        
        MultiArg<string> ideal_bounds_arg("I", "ideal-bound", "Bounds the given objective function by factor times the ideal heuristic value of this objective function.", false,
                                          "objective:factor");

        MultiArg<string> absolute_bounds_arg("i", "absolute-bound", "Bounds the given objective function by the given value.", false, "objective:value");

        cmd.add(epsilon_argument);
        cmd.add(epsilon_all_argument);
        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(ideal_bounds_arg);
        cmd.add(absolute_bounds_arg);
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
        
        calculate_ideal_heuristic(graph,
                                  costs,
                                  dimension,
                                  source,
                                  target,
                                  directed,
                                  distances);
        
        auto ideal_heuristic = [distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        Point bounds(numeric_limits<double>::infinity(), dimension);
        parse_ideal_bounds(ideal_bounds_arg,
                           dimension,
                           ideal_heuristic,
                           source,
                           bounds);

        parse_absolute_bounds(absolute_bounds_arg,
                              dimension,
                              bounds);
        
        auto cost_function = [&costs] (edge e) -> Point& {
            return costs[e];
        };

        ConstrainedReachabilityPreprocessing prepro;
        list<Point> bounds_list;
        for(unsigned i = 0; i < dimension; ++i) {
            if(bounds[i] < numeric_limits<double>::infinity()) {
                Point new_bound(dimension + 1);
                new_bound[i] = 1;
                new_bound[dimension] = -bounds[i];
                bounds_list.push_back(std::move(new_bound));
            }
        }
        prepro.preprocess(graph,
                          cost_function,
                          dimension,
                          source,
                          target,
                          bounds_list);

        
        Point epsilon(epsilon_all, dimension);
        if(exact_argument.isSet()) {
            epsilon[exact_objective - 1] = 0.0;
        }
        
        parse_epsilon(epsilon_argument,
                      dimension,
                      epsilon);
        
        LCApprox solver;
        
        solver.set_bounds(bounds);
        solver.set_heuristic(ideal_heuristic);
        
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

void EpLCApproxModule::parse_ideal_bounds(const MultiArg<string>& argument,
                                         unsigned dimension,
                                         function<double(node, unsigned)> heuristic,
                                         const node source,
                                         Point& bounds) {
    
    auto bounds_it = argument.begin();
    while(bounds_it != argument.end()) {
        vector<string> tokens;
        mco::tokenize(*bounds_it, tokens, ":");
        
        if(tokens.size() != 2) {
            cout << "Error" << endl;
            return;
        }
        
        unsigned objective_function = stoul(tokens[0]);
        double factor = stod(tokens[1]);
        
        if(objective_function > dimension) {
            cout << "Error" << endl;
            return;
        }
        
        if(objective_function == 0) {
            cout << "Error" << endl;
            return;
        }
        
        bounds[objective_function - 1] = factor * heuristic(source, objective_function - 1);
        
        bounds_it++;
    }
    
    
}

void EpLCApproxModule::parse_absolute_bounds(const MultiArg<string>& argument,
                                          unsigned dimension,
                                          Point& bounds) {

    auto bounds_it = argument.begin();
    while(bounds_it != argument.end()) {
        vector<string> tokens;
        mco::tokenize(*bounds_it, tokens, ":");

        if(tokens.size() != 2) {
            cout << "Error" << endl;
            return;
        }

        unsigned objective_function = stoul(tokens[0]);
        double value = stod(tokens[1]);

        if(objective_function > dimension) {
            cout << "Error" << endl;
            return;
        }

        if(objective_function == 0) {
            cout << "Error" << endl;
            return;
        }

        bounds[objective_function - 1] = value;

        bounds_it++;
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

void EpLCApproxModule::calculate_ideal_heuristic(const Graph& graph,
                                                 const EdgeArray<Point>& costs,
                                                 unsigned dimension,
                                                 const node source,
                                                 const node target,
                                                 bool directed,
                                                 vector<NodeArray<double>>& distances) {

    Dijkstra<double> sssp_solver;

    NodeArray<edge> predecessor(graph);

    for(unsigned i = 0; i < dimension; ++i) {
        auto length = [&costs, i] (edge e) {
            return costs[e][i];
        };

        sssp_solver.singleSourceShortestPaths(graph,
                                              length,
                                              target,
                                              predecessor,
                                              distances[i],
                                              directed ? DijkstraModes::Backward :
                                              DijkstraModes::Undirected);
    }

}


const list<pair<const list<edge>, const Point>>& EpLCApproxModule::solutions() {
    return solutions_;
}

string EpLCApproxModule::statistics() {
    string stats("");
    return stats;
}
