//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "ep_ideal_point.h"

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

using mco::EPDualBensonSolver;
using mco::TemporaryGraphParser;
using mco::Point;
using mco::EpSolverMartins;
using mco::Dijkstra;
using mco::DijkstraModes;
using mco::ConstrainedReachabilityPreprocessing;

void EpIdealModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Computes the ideal point the efficient path problem.", ' ', "0.1");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);
        

        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        bool is_directed = is_directed_arg.getValue();

        Graph graph;
        EdgeArray<Point> costs(graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, graph, costs, dimension, source, target);

        Point ideal_point(dimension);

        ogdf::NodeArray<double> distances(graph);
        ogdf::NodeArray<edge> predecessors(graph);

        for(unsigned i = 0; i < dimension; ++i) {
            mco::Dijkstra<double, mco::PairComparator<double, std::less<double>>> sssp;

            auto cost_function = [&costs, i] (edge e) {
                return costs[e][i];
            };

            sssp.singleSourceShortestPaths(graph,
                                           cost_function,
                                           source,
                                           predecessors,
                                           distances,
                                           is_directed ? DijkstraModes::Forward : DijkstraModes::Undirected);

            ideal_point[i] = distances[target];
        }

        cout << ideal_point << endl;
        

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EpIdealModule::solutions() {
    return solutions_;
}

string EpIdealModule::statistics() {
    return string("");
}