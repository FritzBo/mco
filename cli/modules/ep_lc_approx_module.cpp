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

void EpLCApproxModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Label Correcting Approximation for the Efficient Path Problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "The approximation factor to use.", true, 1, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);
        
        MultiArg<string> ideal_bounds_arg("I", "ideal-bound", "Bounds the given objective function by factor times the ideal heuristic value of this objective function. Implies -H.", false,
                                          "objective:factor");
        
        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(ideal_bounds_arg);
        
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
        
        list<Point> bounds;
        parse_ideal_bounds(ideal_bounds_arg,
                           dimension,
                           ideal_heuristic,
                           source,
                           bounds);
        
        auto cost_function = [&costs] (edge e) -> Point& {
            return costs[e];
        };
        
        ConstrainedReachabilityPreprocessing prepro;
        
        prepro.preprocess(graph,
                          cost_function,
                          dimension,
                          source,
                          target,
                          bounds);
        
        
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

void EpLCApproxModule::parse_ideal_bounds(const MultiArg<string>& argument,
                                         unsigned dimension,
                                         function<double(node, unsigned)> heuristic,
                                         const node source,
                                         list<Point>& bounds) {
    
    auto bounds_it = argument.begin();
    while(bounds_it != argument.end()) {
        Point bound(0.0, dimension + 1);
        
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
        
        bound[objective_function - 1] = 1;
        bound[dimension] = -factor * heuristic(source, objective_function - 1);
        
        bounds.push_back(std::move(bound));
        
        bounds_it++;
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
