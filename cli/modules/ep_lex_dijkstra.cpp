//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "ep_lex_dijkstra.h"

#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::list;
using std::pair;
using std::make_pair;
using std::function;
using std::vector;
using std::cout;
using std::endl;

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/Thread.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;
using ogdf::Thread;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ArgException;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;
using TCLAP::MultiArg;

#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/martins/martins.h>
#include <mco/ep/preprocessing/constrained_reach.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>
#include <mco/cli/parse_util.h>

using mco::EPDualBensonSolver;
using mco::TemporaryGraphParser;
using mco::Point;
using mco::EpSolverMartins;
using mco::LexDijkstra;
using mco::DijkstraModes;
using mco::ConstrainedReachabilityPreprocessing;

void EpLDModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Computes a lexicographic minimal path with optional weights.", ' ', "0.1");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);

        ValueArg<string> weight_arg("W", "weight", "Adding a weighted-sum as first objective", false, "", "weight");
        

        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(weight_arg);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        bool is_directed = is_directed_arg.getValue();

        Graph graph;
        EdgeArray<Point> costs(graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, graph, costs, dimension, source, target);

        ogdf::NodeArray<Point*> distances(graph);
        ogdf::NodeArray<edge> predecessors(graph);

        LexDijkstra sssp;

        if(weight_arg.isSet()) {

            Point weight = parse_weight(weight_arg, dimension);

            EdgeArray<Point> weighted_cost(graph);

            for(auto e : graph.edges) {
                weighted_cost[e] = Point(dimension + 1);

                weighted_cost[e][0] = weight * costs[e];

                for(unsigned i = 1; i <= dimension; ++i) {
                    weighted_cost[e][i] = costs[e][i - 1];
                }
            }

            for(auto n : graph.nodes) {
                distances[n] = new Point(dimension + 1);
            }

            auto cost_function = [&weighted_cost] (edge e) {
                return &weighted_cost[e];
            };

            sssp.singleSourceShortestPaths(graph,
                                           cost_function,
                                           source,
                                           distances,
                                           predecessors,
                                           is_directed ? DijkstraModes::Forward : DijkstraModes::Undirected);



        } else {
            auto cost_function = [&costs] (edge e) {
                return &costs[e];
            };

            for(auto n : graph.nodes) {
                distances[n] = new Point(dimension);
            }

            sssp.singleSourceShortestPaths(graph,
                                           cost_function,
                                           source,
                                           distances,
                                           predecessors,
                                           is_directed ? DijkstraModes::Forward : DijkstraModes::Undirected);

        }

        cout << *distances[target] << endl;

        list<node> path;
        auto n = target;
        path.push_front(n);

        while(n != source) {

            n = predecessors[n]->opposite(n);
            path.push_front(n);
        }

        for(auto n : path) {
            cout << n << ", ";
        }

        for(auto n : graph.nodes) {
            delete distances[n];
        }

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EpLDModule::solutions() {
    return solutions_;
}

string EpLDModule::statistics() {
    return string("");
}

Point EpLDModule::parse_weight(ValueArg<string>& argument,
                               unsigned dimension) {

    vector<string> tokens;
    mco::tokenize(argument.getValue(), tokens, ":");

    if(tokens.size() != dimension) {
        cout << "Error: Wrong number of weights" << endl;
        return Point();
    }

    Point r(dimension);

    for(unsigned i = 0; i < dimension; ++i) {
        r[i] = stod(tokens[i]);
    }

    return r;
}