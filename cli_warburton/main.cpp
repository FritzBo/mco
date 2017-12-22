//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 04.08.14.
//
//

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <chrono>

using std::list;
using std::pair;
using std::vector;
using std::string;
using std::numeric_limits;
using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::ValueArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;
using TCLAP::MultiArg;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;

#include <mco/basic/point.h>

#include <mco/benchmarks/temporary_graphs_parser.h>

#include <mco/cli/parse_util.h>
#include <mco/cli/output_formatter.h>

#include <mco/ep/basic/dijkstra.h>
#include <mco/ep/basic/instance_scalarizer.h>
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
        
        SwitchArg directed_arg("d",
                               "directed",
                               "Switch if the input should be interpreted as "
                               "being directed.",
                               false);
        
        SwitchArg test_only_arg("t",
                                "test-only",
                                "Perform only a test of how many subproblemes "
                                "will be created.",
                                false);
        
        MultiArg<string> ideal_bounds_arg("i",
                                          "ideal-bound",
                                          "Bounds the given objective function "
                                          "by a factor times the ideal value "
                                          "of this objective function",
                                          false,
                                          "objective:factor");

        SwitchArg print_front_arg("f",
        "print-fronrt",
        "Print the whole nondominated set instead of statistics only.", false);
            
        cmd_line.add(filename_arg);
        cmd_line.add(epsilon_arg);
        cmd_line.add(do_preprocessing_arg);
        cmd_line.add(directed_arg);
        cmd_line.add(test_only_arg);
        cmd_line.add(ideal_bounds_arg);
        cmd_line.add(print_front_arg);
        
        CliOutputFormatter<list<ogdf::edge>> output_formatter(&cmd_line);
        
        output_formatter.initialize_cmd_line();
        
        cmd_line.parse(argc, argv);
            
        string filename         = filename_arg.getValue();
        double epsilon          = epsilon_arg.getValue();
        bool do_preprocessing   = do_preprocessing_arg.getValue();
        bool directed           = directed_arg.getValue();
        bool test_only          = test_only_arg.getValue();
        bool print_front        = print_front_arg.getValue();
            
        Graph graph;
        EdgeArray<Point> raw_costs(graph);
        unsigned dimension;
        ogdf::node source;
        ogdf::node target;
            
        TemporaryGraphParser parser;
        parser.getGraph(filename, graph, raw_costs, dimension, source, target);

        EdgeArray<Point> costs(graph, Point(dimension));
        Point factor(100.0, dimension);
        InstanceScalarizer::scaleround_instance(graph,
                                                raw_costs,
                                                dimension,
                                                factor,
                                                costs);
        
        auto cost_function = [&costs] (ogdf::edge e) -> Point& {
            return costs[e];
        };
        
        auto cost_function2 = [&costs] (ogdf::edge e) {
            return &costs[e];
        };
        
        Point ideal_bound(numeric_limits<double>::infinity(), dimension);
        
        vector<NodeArray<double>> distances(dimension, graph);
        NodeArray<ogdf::edge> predecessor(graph);
        
        for(unsigned i = 0; i < dimension; ++i) {
            Dijkstra<double, PairComparator<double, std::less<double>>> sssp_solver;
            
            auto mono_cost = [&costs, i] (ogdf::edge e) {
                return costs[e][i];
            };
            
            sssp_solver.singleSourceShortestPaths(graph,
                                                  mono_cost,
                                                  target,
                                                  predecessor,
                                                  distances[i],
                                                  directed ? DijkstraModes<>::Backward :
                                                  DijkstraModes<>::Undirected);
        }
        
        auto heuristic = [&distances] (ogdf::node n, unsigned objective) {
            return distances[objective][n];
        };
        
        parse_ideal_bounds(ideal_bounds_arg,
                           dimension,
                           heuristic,
                           source,
                           ideal_bound);
        
        if(do_preprocessing) {
            ConstrainedReachabilityPreprocessing prepro;
            list<Point> linear_bounds;
            for(unsigned i = 0; i < dimension; ++i) {
                if(ideal_bound[i] != numeric_limits<double>::infinity()) {
                    Point linear_bound(dimension + 1);
                    linear_bound[i] = 1;
                    linear_bound[dimension] = -ideal_bound[i];
                    
                    linear_bounds.push_back(std::move(linear_bound));
                }
            }
            prepro.preprocess(graph, cost_function2, dimension, source, target, linear_bounds);
        }

        steady_clock::time_point start = steady_clock::now();
            
        EpSolverWarburtonApprox solver;
            
        solver.Solve(graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     Point(epsilon, dimension),
                     ideal_bound,
                     test_only,
                     directed);

        steady_clock::time_point end = steady_clock::now();
        duration<double> computation_span
                                     = duration_cast<duration<double>>(end - start);

        if(not print_front)
        {
            cout << filename << ", " << graph.numberOfNodes() << ", " << graph.numberOfEdges()
                 << ", " << dimension << ", " << epsilon << ", " << computation_span.count() << ", "
                 << solver.solutions().size() << endl;
        }
        else
        {
            for(auto solution : solver.solutions())
            {
                for(unsigned i = 0; i < solution.second.dimension(); ++i)
                {
                    cout << solution.second[i] << ", ";
                }
                cout << endl;
            }
        }

        
    } catch(TCLAP::ArgException& e) {
        cout << e.typeDescription() << endl;
    }
    
    return 0;
}