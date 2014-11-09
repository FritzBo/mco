//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "ep_benson_module.h"

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

#include <mco/ep/basic/instance_scalarizer.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>

using mco::InstanceScalarizer;
using mco::EPDualBensonSolver;
using mco::OnlineVertexEnumeratorCDD;
using mco::TemporaryGraphParser;
using mco::Point;

void EpBensonModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Dual Benson to find the Pareto-frontier of the efficient path problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 1E-8, "epsilon");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);

        SwitchArg use_cdd_arg("", "cdd", "Using CDD Library for the vertex enumeration", false);

        SwitchArg use_gl_ove_arg("", "gl-ove", "Using the graphless online vertex enumerator", false);
        
        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(use_cdd_arg);
        cmd.add(use_gl_ove_arg);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();
        bool directed = is_directed_arg.getValue();
        bool use_cdd = use_cdd_arg.getValue();
        bool use_gl_ove = use_gl_ove_arg.getValue();
        
        Graph graph;
        EdgeArray<Point> raw_costs(graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, graph, raw_costs, dimension, source, target);

        Point factor(numeric_limits<double>::max(), dimension);
        for(auto e : graph.edges) {
            for(unsigned i = 0; i < dimension; ++i) {
                factor[i] = min(factor[i], 1 / raw_costs[e][i]);
            }
        }

        EdgeArray<Point> costs(graph, Point(dimension));
        InstanceScalarizer::scale_instance(graph,
                                           raw_costs,
                                           dimension,
                                           factor,
                                           costs);
        if((dimension <= 4 && !use_cdd) || use_gl_ove) {
            EPDualBensonSolver<> solver(epsilon);
            
            auto cost_function = [costs] (edge e) { return &costs[e]; };
            
            solver.Solve(graph, cost_function, source, target);

            for(auto sol : solver.solutions()) {
                Point point = sol.second;
                for(unsigned i = 0; i < dimension; ++i) {
                    point[i] /= factor[i];
                }
                solutions_.push_back(make_pair(list<edge>(), point));
            }
        } else {
            EPDualBensonSolver<OnlineVertexEnumeratorCDD> solver(epsilon);

            auto cost_function = [costs] (edge e) { return &costs[e]; };

            solver.Solve(graph, cost_function, source, target);

            for(auto sol : solver.solutions()) {
                Point point = sol.second;
                for(unsigned i = 0; i < dimension; ++i) {
                    point[i] /= factor[i];
                }
                solutions_.push_back(make_pair(list<edge>(), point));
            }
        }
        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EpBensonModule::solutions() {
    return solutions_;
}

string EpBensonModule::statistics() {
    string stats("");
    return stats;
}
