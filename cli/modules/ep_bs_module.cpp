//
//  ep_lc_approx_module.cpp
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#include "ep_bs_module.h"

#include <thread>
#include <sstream>

using std::map;
using std::string;
using std::list;
using std::vector;
using std::pair;
using std::function;
using std::thread;
using std::make_pair;
using std::numeric_limits;
using std::cout;
using std::endl;
using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

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
#include <mco/ep/brum_shier/ep_solver_bs.h>
#include <mco/ep/brum_shier/ep_solver_bs_bi.h>
#include <mco/ep/brum_shier/ep_solver_bs_td.h>
#include <mco/ep/ep_lc_ls.h>
#include <mco/ep/dual_benson/ep_dual_benson.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/basic/point.h>
#include <mco/cli/parse_util.h>

using mco::EpSolverBS;
using mco::EpSolverBSBi;
using mco::EpSolverBSTd;
using mco::EpLcLs;
using mco::TemporaryGraphParser;
using mco::Point;
using mco::Dijkstra;
using mco::DijkstraModes;
using mco::PairComparator;
using mco::ConstrainedReachabilityPreprocessing;
using mco::InstanceScalarizer;
using mco::EPDualBensonSolver;

using mco::ForwardStar;
using mco::node;
using mco::edge;
using mco::FSEdgeArray;
using mco::FSNodeArray;
using mco::ForwardStarFileReader;

// Todo: Delete
#include <chrono>
#include <thread>

void EpBsModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Label Correcting Approximation for the Efficient Path Problem.", ' ', "0.1");

        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);
        
        MultiArg<string> ideal_bounds_arg("I", "ideal-bound", "Bounds the given objective function by factor times the ideal heuristic value of this objective function.", false,
                                          "objective:factor");

        MultiArg<string> absolute_bounds_arg("i", "absolute-bound", "Bounds the given objective function by the given value.", false, "objective:value");

        SwitchArg do_first_phase_argument("f", "first-phase", "Using dual benson in a first phase", false);

        SwitchArg label_select_arg("l", "label-select", "Using a label-selection strategy instead of node selection", false);

        SwitchArg tree_deletion_arg("t", "tree-deletion", "Enabling tree deletion heuristic", false);

        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(ideal_bounds_arg);
        cmd.add(absolute_bounds_arg);
        cmd.add(do_first_phase_argument);
        cmd.add(label_select_arg);
        cmd.add(tree_deletion_arg);

        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        bool directed = is_directed_arg.getValue();
        bool do_first_phase = do_first_phase_argument.getValue();
        bool label_select = label_select_arg.getValue();
        bool tree_deletion = tree_deletion_arg.getValue();

        ForwardStar graph;
        FSNodeArray<int> extern_ids(graph);
        FSEdgeArray<Point> raw_costs(graph);
        unsigned dimension;
        node source, target;
        
        ForwardStarFileReader file_reader;
        file_reader.read(file_name,
                         graph,
                         extern_ids,
                         raw_costs,
                         dimension,
                         source,
                         target,
                         directed);

        num_nodes_ = graph.numberOfNodes();
        num_edges_ = graph.numberOfEdges();
        num_objectives_ = dimension;

//        FSEdgeArray<Point> costs(graph, Point(dimension));
//        Point factor(100.0, dimension);
//        InstanceScalarizer::scaleround_instance(graph,
//                                                raw_costs,
//                                                dimension,
//                                                factor,
//                                                costs);


//        vector<NodeArray<double>> distances(dimension, graph);
//        
//        calculate_ideal_heuristic(graph,
//                                  costs,
//                                  dimension,
//                                  source,
//                                  target,
//                                  directed,
//                                  distances);
//        
//        auto ideal_heuristic = [distances] (node n, unsigned objective) {
//            return distances[objective][n];
//        };

//        Point bounds(numeric_limits<double>::infinity(), dimension);
//        parse_ideal_bounds(ideal_bounds_arg,
//                           dimension,
//                           ideal_heuristic,
//                           source,
//                           bounds);

//        parse_absolute_bounds(absolute_bounds_arg,
//                              dimension,
//                              bounds);

        auto cost_function = [&raw_costs] (edge e) -> Point* {
            return &raw_costs[e];
        };

//        ConstrainedReachabilityPreprocessing prepro;
//        list<Point> bounds_list;
//        for(unsigned i = 0; i < dimension; ++i) {
//            if(bounds[i] < numeric_limits<double>::infinity()) {
//                Point new_bound(dimension + 1);
//                new_bound[i] = 1;
//                new_bound[dimension] = -bounds[i];
//                bounds_list.push_back(std::move(new_bound));
//            }
//        }
//        prepro.preprocess(graph,
//                          cost_function,
//                          dimension,
//                          source,
//                          target,
//                          bounds_list);


//        thread *first_phase_thread = nullptr;

//        if(do_first_phase) {
//            cout << "starting first phase..." << endl;
//
//            scaled_costs = new EdgeArray<Point>(graph, Point(dimension));
//            Point factor(numeric_limits<double>::max(), dimension);
//            for(auto e : graph.edges) {
//                for(unsigned i = 0; i < dimension; ++i) {
//                    factor[i] = min(factor[i], 1 / raw_costs[e][i]);
//                }
//            }
//
//            InstanceScalarizer::scale_instance(graph,
//                                               raw_costs,
//                                               dimension,
//                                               factor,
//                                               *scaled_costs);
//
//            auto p_cost_function = [scaled_costs] (edge e) {
//                return &scaled_costs->operator[](e);
//            };

//            auto callback = [&solver, &graph, dimension] (NodeArray<Point *>& distances,
//                                                          NodeArray<edge>& predecessors) {
//
//                NodeArray<Point *> new_distances(graph);
//
//                for(auto n : distances.graphOf()->nodes) {
//                    Point* p = new Point(dimension);
//                    std::copy(distances[n]->cbegin() + 1, distances[n]->cend(),
//                              p->begin());
//                    new_distances[n] = p;
//                }
//
//                solver.add_pending_labels(make_pair(new_distances, predecessors));
//            };

//            first_phase_thread = new thread(EpLCApproxModule::first_phase,
//                                            &graph,
//                                            p_cost_function,
//                                            dimension,
//                                            source,
//                                            target);
//
//        }

//        EdgeArray<Point>* scaled_costs = nullptr;

//        if(dimension == 2)
//        {
//
////            solver.set_bounds(bounds);
////            solver.set_heuristic(ideal_heuristic);
//            EpSolverBSBi solver;
//
//            steady_clock::time_point start = steady_clock::now();
//
//            solver.Solve(graph,
//                         cost_function,
//                         dimension,
//                         source,
//                         target);
//
//            steady_clock::time_point end = steady_clock::now();
//            duration<double> computation_span
//            = duration_cast<duration<double>>(end - start);
//            solution_time_ = computation_span.count();
//
//            solutions_.insert(solutions_.begin(),
//                              solver.solutions().cbegin(),
//                              solver.solutions().cend());
//
//            label_compares_ = solver.label_compares();
//            deleted_tree_labels_ = solver.deleted_tree_labels();
//            recursive_deletions_ = solver.recursive_deletions();
//            arc_pushes_ = solver.arc_pushes();
//            touched_recursively_deleted_label_ = solver.touched_recursively_deleted_label();
//            deleted_labels_ = solver.deleted_labels();
//            method_ = "ns-bi";
//        }
//        else
        if(!label_select)
        {
            if(tree_deletion)
            {
                EpSolverBSTd solver;

                steady_clock::time_point start = steady_clock::now();

                solver.Solve(graph,
                             cost_function,
                             dimension,
                             source,
                             target);

                steady_clock::time_point end = steady_clock::now();
                duration<double> computation_span
                = duration_cast<duration<double>>(end - start);
                solution_time_ = computation_span.count();

                solutions_.insert(solutions_.begin(),
                                  solver.solutions().cbegin(),
                                  solver.solutions().cend());

                label_compares_ = solver.label_compares();
                deleted_tree_labels_ = solver.deleted_tree_labels();
                recursive_deletions_ = solver.recursive_deletions();
                arc_pushes_ = solver.arc_pushes();
                touched_recursively_deleted_label_ = solver.touched_recursively_deleted_label();
                deleted_labels_ = solver.deleted_labels();

                method_ = "ns-td";
            }
            else
            {
                EpSolverBS solver;
                
//            solver.set_bounds(bounds);
//            solver.set_heuristic(ideal_heuristic);

                steady_clock::time_point start = steady_clock::now();

                solver.Solve(graph,
                             cost_function,
                             dimension,
                             source,
                             target);

                steady_clock::time_point end = steady_clock::now();
                duration<double> computation_span
                = duration_cast<duration<double>>(end - start);
                solution_time_ = computation_span.count();

                solutions_.insert(solutions_.begin(),
                                  solver.solutions().cbegin(),
                                  solver.solutions().cend());

                label_compares_ = solver.label_compares();
                deleted_tree_labels_ = solver.deleted_tree_labels();
                recursive_deletions_ = solver.recursive_deletions();
                arc_pushes_ = solver.arc_pushes();
                touched_recursively_deleted_label_ = solver.touched_recursively_deleted_label();
                deleted_labels_ = solver.deleted_labels();
                method_ = "ns";
            }
        }
        else
        {
            EpLcLs solver;


//            solver.set_bounds(bounds);
//            solver.set_heuristic(ideal_heuristic);

            solver.Solve(graph,
                         cost_function,
                         dimension,
                         source,
                         target);

            solutions_.insert(solutions_.begin(),
                              solver.solutions().cbegin(),
                              solver.solutions().cend());

            label_compares_ = solver.label_compares();
            deleted_tree_labels_ = solver.deleted_tree_labels();
            recursive_deletions_ = solver.recursive_deletions();
            arc_pushes_ = solver.arc_pushes();
            touched_recursively_deleted_label_ = solver.touched_recursively_deleted_label();
            deleted_labels_ = solver.deleted_labels();
            method_ = "ls";
        }

//        delete scaled_costs;

    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

//void EpBsModule::parse_ideal_bounds(const MultiArg<string>& argument,
//                                         unsigned dimension,
//                                         function<double(node, unsigned)> heuristic,
//                                         const node source,
//                                         Point& bounds) {
//    
//    auto bounds_it = argument.begin();
//    while(bounds_it != argument.end()) {
//        vector<string> tokens;
//        mco::tokenize(*bounds_it, tokens, ":");
//        
//        if(tokens.size() != 2) {
//            cout << "Error" << endl;
//            return;
//        }
//        
//        unsigned objective_function = stoul(tokens[0]);
//        double factor = stod(tokens[1]);
//        
//        if(objective_function > dimension) {
//            cout << "Error" << endl;
//            return;
//        }
//        
//        if(objective_function == 0) {
//            cout << "Error" << endl;
//            return;
//        }
//        
//        bounds[objective_function - 1] = factor * heuristic(source, objective_function - 1);
//        
//        bounds_it++;
//    }
//    
//    
//}

void EpBsModule::parse_absolute_bounds(const MultiArg<string>& argument,
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

//void EpBsModule::calculate_ideal_heuristic(const Graph& graph,
//                                                 const EdgeArray<Point>& costs,
//                                                 unsigned dimension,
//                                                 const node source,
//                                                 const node target,
//                                                 bool directed,
//                                                 vector<NodeArray<double>>& distances) {
//
//    Dijkstra<double, PairComparator<double, std::less<double>>> sssp_solver;
//
//    NodeArray<edge> predecessor(graph);
//
//    for(unsigned i = 0; i < dimension; ++i) {
//        auto length = [&costs, i] (edge e) {
//            return costs[e][i];
//        };
//
//        sssp_solver.singleSourceShortestPaths(graph,
//                                              length,
//                                              target,
//                                              predecessor,
//                                              distances[i],
//                                              directed ? DijkstraModes::Backward :
//                                              DijkstraModes::Undirected);
//    }
//
//}

//void EpBsModule::first_phase(const Graph* graph,
//                                  function<Point*(edge)> cost_function,
//                                  unsigned dimension,
//                                  const node source,
//                                  const node target,
//                                  std::function<void(ogdf::NodeArray<Point*>&, ogdf::NodeArray<ogdf::edge>&)> callback) {
//
//    OGDF_ALLOCATOR::initThread();
//    //    for(double epsilon = 1; epsilon > 1E-6; epsilon /= 10) {
//    EPDualBensonSolver<> weighted_solver(1E-8);
//
//    weighted_solver.Solve(*graph, cost_function, source, target, callback);
//    //    }
//    OGDF_ALLOCATOR::flushPool();
//
//    cout << "Finished first phase" << endl;
//}

const list<pair<const list<mco::node>, const Point>>& EpBsModule::
solutions() {
    return solutions_;
}

string EpBsModule::statistics() {
    std::stringstream sstream;
    sstream << num_nodes_;
    sstream << ", ";
    sstream << num_edges_;
    sstream << ", ";
    sstream << num_objectives_;
    sstream << ", ";
    sstream << method_;
    sstream << ", ";
    sstream << label_compares_;
    sstream << ", ";
    sstream << deleted_labels_;
    sstream << ", ";
    sstream << arc_pushes_;
    sstream << ", ";
    sstream << deleted_tree_labels_;
    sstream << ", ";
    sstream << recursive_deletions_;
    sstream << ", ";
    sstream << touched_recursively_deleted_label_;
    sstream << ", ";
    sstream << solution_time_;

    return sstream.str();
}
