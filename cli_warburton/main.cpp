//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 04.08.14.
//
//

#include <list>

using std::list;
using std::pair;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::node;

#include <mco/basic/point.h>

#include <mco/benchmarks/temporary_graphs_parser.h>

#include <mco/cli/output_formatter.h>

#include <mco/ep/preprocessing/constrained_reach.h>
#include <mco/ep/warburton/ep_solver_warburton_approx.h>

using namespace mco;

int main(int argc, char** argv) {
    
    try {
    
        CmdLine cmd_line("Implementation of Warburton' approximation algorithm "
                         "for the multiobjective shortest s-t-path problem",
                         ' ',
                         "0.1");
            
        ValueArg<double> epsilon_arg("e",
                                     "tolerance",
                                     "Floating point value to control the "
                                     "tolerance of the algorithm",
                                     true,
                                     0,
                                     "float");
            
        SwitchArg do_preprocessing_arg("p",
                                       "preprocessing",
                                       "Conduct a first preprocessing step",
                                       false);
            
        UnlabeledValueArg<string> filename_arg("filename",
                                               "Name of the instance file",
                                               true,
                                               "",
                                               "filename");
            
        cmd_line.add(filename_arg);
        cmd_line.add(epsilon_arg);
        cmd_line.add(do_preprocessing_arg);
        
        CliOutputFormatter<list<edge>> output_formatter(&cmd_line);
        
        output_formatter.initialize_cmd_line();
        
        cmd_line.parse(argc, argv);
            
        string filename         = filename_arg.getValue();
        double epsilon          = epsilon_arg.getValue();
        bool do_preprocessing   = do_preprocessing_arg.getValue();
            
        Graph graph;
        EdgeArray<Point> costs(graph);
        unsigned dimension;
        node source;
        node target;
            
        TemporaryGraphParser parser;
        parser.getGraph(filename, graph, costs, dimension, source, target);
        
        auto cost_function = [&costs] (edge e) -> Point& {
            return costs[e];
        };
        
        auto cost_function2 = [&costs] (edge e) {
            return &costs[e];
        };

        if(do_preprocessing) {
            ConstrainedReachabilityPreprocessing prepro;
            Point length_bound({0, 0, 0, 0, 0, 1, -34000});
            Point bundling_bound({1, 0, 0, 0, 0, 0, -26000});
            list<Point> linear_bounds({length_bound, bundling_bound});
            prepro.preprocess(graph, cost_function2, dimension, source, target, linear_bounds);
        }

            
        EpSolverWarburtonApprox solver;
            
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     Point(epsilon, dimension));
        
        output_formatter.print_output(cout,
                                      solver.solutions().begin(),
                                      solver.solutions().end());
        
    } catch(TCLAP::ArgException& e) {
        cout << e.typeDescription() << endl;
    }
    
    return 0;
}