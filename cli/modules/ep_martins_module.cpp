//
//  benson_module.cpp
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#include "ep_martins_module.h"

#include <map>
#include <string>
#include <vector>
#include <thread>

using std::map;
using std::string;
using std::list;
using std::pair;
using std::make_pair;
using std::function;
using std::vector;
using std::thread;

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

void EpMartinsModule::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Martins Algorithm to find the Pareto-frontier of the efficient path problem.", ' ', "0.1");
        
        ValueArg<double> epsilon_argument("e", "epsilon", "Epsilon to be used in floating point calculations.", false, 0, "epsilon");
        
//        ValueArg<string> bundling_arg("B", "bundling", "", false, "", "objective:...:objective:reference_objective:factor");
        
        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");
        
        SwitchArg use_heuristic_switch("H", "use-heuristic", "Use the ideal point heuristic to prune paths", false);
        
        SwitchArg do_first_phase_arg("f", "first-phase", "Conducts a first phase of finding extremal points of the Pareto-frontier. Only useful if computation takes very long.", false);
        
        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);
        
        MultiArg<string> ideal_bounds_arg("I", "ideal-bound", "Bounds the given objective function by factor times the ideal heuristic value of this objective function. Implies -H.", false,
                                     "objective:factor");
        
//        MultiArg<string> fractional_bounds_arg("F", "frational-bound", "Bounds the quotients of the given objective functions by a factor. Implies -H.", false, "bound_objective:reference_objective:factor");
        
        cmd.add(epsilon_argument);
        cmd.add(file_name_argument);
        cmd.add(use_heuristic_switch);
        cmd.add(ideal_bounds_arg);
//        cmd.add(fractional_bounds_arg);
//        cmd.add(bundling_arg);
        cmd.add(is_directed_arg);
        cmd.add(do_first_phase_arg);
        
        cmd.parse(argc, argv);
        
        string file_name = file_name_argument.getValue();
        double epsilon = epsilon_argument.getValue();
        bool use_heuristic = use_heuristic_switch.getValue();
        bool is_directed = is_directed_arg.getValue();
        bool do_first_phase = do_first_phase_arg.getValue();
        
        if(ideal_bounds_arg.end() - ideal_bounds_arg.begin() > 0) {
//           fractional_bounds_arg.end() - fractional_bounds_arg.begin() > 0 ||
//           bundling_arg.isSet()) {
           
            use_heuristic = true;
        }
        
        Graph* graph = new Graph;
        EdgeArray<Point> costs(*graph);
        unsigned dimension;
        node source, target;
        
        TemporaryGraphParser parser;
        
        parser.getGraph(file_name, *graph, costs, dimension, source, target);
        
        vector<NodeArray<double>> distances(dimension, *graph);
        
        if(use_heuristic) {
#ifndef NDEBUG
            cout << "computing heuristic..." <<endl;
#endif
            calculate_ideal_heuristic(*graph,
                                      costs,
                                      dimension,
                                      source,
                                      target,
                                      is_directed,
                                      distances);
            
        } else {
            for(unsigned i = 0; i < dimension; ++i) {
                for(auto n : graph->nodes) {
                    distances[i][n] = 0;
                }
            }
        }
        
        auto ideal_heuristic = [distances] (node n, unsigned objective) {
            return distances[objective][n];
        };
        
        list<Point> bounds;
        parse_ideal_bounds(ideal_bounds_arg,
                           dimension,
                           ideal_heuristic,
                           source,
                           bounds);
        
//        parse_fractional_bounds(fractional_bounds_arg,
//                                dimension,
//                                bounds);
        
//        if(bundling_arg.isSet()) {
//            parse_bundling_bound(bundling_arg,
//                                 dimension,
//                                 ideal_heuristic,
//                                 bounds,
//                                 bundling_bound);
//        }
        
        auto cost_function = [&costs] (edge e) { return &costs[e]; };
        
        ConstrainedReachabilityPreprocessing prep;

#ifndef NDEBUG
        cout << "starting preprocessing..." << endl;
#endif
        prep.preprocess(*graph,
                        cost_function,
                        dimension,
                        source,
                        target,
                        bounds,
                        is_directed);
        
        
        EpSolverMartins solver(epsilon);
        
        
        thread *first_phase_thread = nullptr;
        if(do_first_phase) {
#ifndef NDEBUG
            cout << "starting first phase..." << endl;
#endif
            
            auto callback = [&solver, &graph, dimension] (NodeArray<Point *>& distances,
                                                          NodeArray<edge>& predecessors) {
                
                NodeArray<Point *> new_distances(*graph);
                
                for(auto n : distances.graphOf()->nodes) {
                    Point* p = new Point(dimension);
                    std::copy(distances[n]->cbegin() + 1, distances[n]->cend(),
                              p->begin());
                    new_distances[n] = p;
                }
                
                solver.add_pending_labels(make_pair(new_distances, predecessors));
            };

            first_phase_thread = new thread(EpMartinsModule::first_phase,
                                            graph,
                                            cost_function,
                                            dimension,
                                            source,
                                            target,
                                            epsilon,
                                            callback);
            
        }

#ifndef NDEBUG
        cout << "Starting Martins..." << endl;
#endif
        
        solver.Solve(*graph,
                     cost_function,
                     dimension,
                     source,
                     target,
                     Point(numeric_limits<double>::infinity(), dimension),
                     std::list<std::pair<ogdf::NodeArray<Point*>, ogdf::NodeArray<ogdf::edge>>>(),
                     ideal_heuristic,
                     bounds,
                     is_directed);

#ifndef NDEBUG
        cout << "Done." << endl;
#endif
        
        if(first_phase_thread != nullptr) {
            first_phase_thread->join();
            delete first_phase_thread;
        }
        
        delete graph;

        solutions_.insert(solutions_.begin(),
                          solver.solutions().cbegin(),
                          solver.solutions().cend());
        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

const list<pair<const list<edge>, const Point>>& EpMartinsModule::solutions() {
    return solutions_;
}

string EpMartinsModule::statistics() {
    return string("");
}

void EpMartinsModule::parse_ideal_bounds(const MultiArg<string>& argument,
                                         unsigned dimension,
                                         function<double(node, unsigned)> heuristic,
                                         const node source,
                                         list<Point>& bounds) {
    
    auto bounds_it = argument.begin();
    while(bounds_it != argument.end()) {
        Point bound(0.0, dimension + 1);
        
        vector<string> tokens;
        tokenize(*bounds_it, tokens, ":");
        
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

void EpMartinsModule::parse_fractional_bounds(const MultiArg<string>& argument,
                                              unsigned dimension,
                                              Point& bounds) {
    
    auto bounds_it = argument.begin();
    while(bounds_it != argument.end()) {
        vector<string> tokens;
        tokenize(*bounds_it, tokens, ":");
        
        if(tokens.size() != 3) {
            cout << "Error" << endl;
            return;
        }
        
        unsigned bound_objective = stoul(tokens[0]);
        unsigned ref_objective = stoul(tokens[1]);
        double factor = stod(tokens[2]);
        
        if(bound_objective > dimension ||
           ref_objective > dimension) {
            
            cout << "Error" << endl;
            return;
        }
        
        if(bound_objective == 0 ||
           ref_objective == 0) {
            
            cout << "Error" << endl;
            return;
        }
        
        if(bounds[ref_objective - 1] == numeric_limits<double>::infinity()) {
            cout << "Error" << endl;
            return;
        }
        
        bounds[bound_objective - 1] = bounds[ref_objective - 1] * factor;
        
        bounds_it++;
    }
    
}

void EpMartinsModule::parse_bundling_bound(ValueArg<string>& bundling_arg,
                                           unsigned dimension,
                                           function<double(node, unsigned)> ideal_heuristic,
                                           const Point& ideal_bound,
                                           Point& bundling_bound) {
    
    vector<string> tokens;
    tokenize(bundling_arg.getValue(), tokens, ":");
    
    if(tokens.size() < 3) {
        cout << "Error: Bounding needs at least 3 parameter." << endl;
        return;
    }
    
    vector<unsigned> bundling_objectives;
    for(unsigned i = 0; i < tokens.size() - 2; ++i) {
        bundling_objectives.push_back(stoul(tokens[i]));
    }
    std::sort(bundling_objectives.begin(),
              bundling_objectives.end());
    
    unsigned reference_objective = stoul(tokens[tokens.size() - 2]);
    double factor = stod(tokens[tokens.size() - 1]);
    
    if(ideal_bound[reference_objective - 1] == numeric_limits<double>::infinity()) {
        cout << "Error: Reference objective must be bounded." << endl;
        return;
    }
    
    auto objective_it = bundling_objectives.begin();
    for(unsigned i = 0; i < dimension; ++i) {
        if(objective_it != bundling_objectives.end() &&
           i == (*objective_it) - 1) {
            
            bundling_bound[i] = 1;
            objective_it++;

        } else {
            bundling_bound[i] = 0;
        }
    }
    
    bundling_bound[dimension] = -(bundling_objectives.size() - factor) * ideal_bound[reference_objective - 1];
}


void EpMartinsModule::first_phase(const Graph* graph,
                                  function<Point*(edge)> cost_function,
                                  unsigned dimension,
                                  const node source,
                                  const node target,
                                  double epsilon,
                                  std::function<void(ogdf::NodeArray<Point*>&, ogdf::NodeArray<ogdf::edge>&)> callback) {
    
    OGDF_ALLOCATOR::initThread();
//    for(double epsilon = 1; epsilon > 1E-6; epsilon /= 10) {
        EPDualBensonSolver<> weighted_solver(1E-8);
            
        weighted_solver.Solve(*graph, cost_function, source, target, callback);
//    }
    OGDF_ALLOCATOR::flushPool();
    
    cout << "Finished first phase" << endl;
}


void EpMartinsModule::calculate_ideal_heuristic(
                    const Graph& graph,
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

/*
 } else if(algorithm.compare("pre-martins") == 0) {
 Dijkstra<double> sssp_solver;
 
 vector<NodeArray<double>> distances(dimension, graph);
 NodeArray<edge> predecessor(graph);
 
 cout << "calculating heuristic..." << endl;
 
 for(unsigned i = 0; i < dimension; ++i) {
 auto length = [&costs, i] (edge e) {
 return costs[e][i];
 };
 
 sssp_solver.singleSourceShortestPaths(graph,
 length,
 target,
 predecessor,
 distances[i],
 DijkstraModes::Undirected);
 }
 
 auto heuristic = [&distances] (node n, unsigned objective) {
 return distances[objective][n];
 };
 
 Point absolute_bound(numeric_limits<double>::infinity(), dimension);
 for(unsigned i = 0; i < dimension; ++i) {
 absolute_bound[i] = numeric_limits<double>::infinity();
 }
 //        bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
 
 EpSolverMartins solver;
 
 cout << "Running Martins algorithm..." << endl;
 
 solver.Solve(graph,
 cost_function,
 dimension,
 source,
 target,
 absolute_bound,
 heuristic,
 list<Point>(),
 false);
 
 cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
 
 } else if(algorithm.compare("fpre-martins") == 0) {
 Dijkstra<double> sssp_solver;
 
 vector<NodeArray<double>> distances(dimension, graph);
 NodeArray<edge> predecessor(graph);
 
 cout << "calculating heuristic..." << endl;
 
 for(unsigned i = 0; i < dimension; ++i) {
 auto length = [&costs, i] (edge e) {
 return costs[e][i];
 };
 
 sssp_solver.singleSourceShortestPaths(graph,
 length,
 target,
 predecessor,
 distances[i],
 DijkstraModes::Undirected);
 }
 
 auto heuristic = [&distances] (node n, unsigned objective) {
 return distances[objective][n];
 };
 
 Point absolute_bound(numeric_limits<double>::infinity(), dimension);
 for(unsigned i = 0; i < dimension - 1; ++i) {
 absolute_bound[i] = numeric_limits<double>::infinity();
 }
 absolute_bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
 
 cout << "Running first phase..." << endl;
 
 list<Point> first_phase_bounds;
 {
 
 EPDualBensonSolver<> weighted_solver;
 
 weighted_solver.Solve(graph, cost_function, source, target);
 
 
 for(auto point : weighted_solver.solutions()) {
 first_phase_bounds.push_back(*point);
 }
 
 }
 
 EpSolverMartins solver;
 
 cout << "Running Martins algorithm..." << endl;
 
 solver.Solve(graph,
 cost_function,
 dimension,
 source,
 target,
 absolute_bound,
 heuristic,
 first_phase_bounds,
 false);
 
 cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
 
 } else if(algorithm.compare("flabel-martins") == 0) {
 Dijkstra<double> sssp_solver;
 
 vector<NodeArray<double>> distances(dimension, graph);
 NodeArray<edge> predecessor(graph);
 
 cout << "Calculating heuristic... ";
 steady_clock::time_point start = steady_clock::now();
 
 for(unsigned i = 0; i < dimension; ++i) {
 auto length = [&costs, i] (edge e) {
 return costs[e][i];
 };
 
 sssp_solver.singleSourceShortestPaths(graph,
 length,
 target,
 predecessor,
 distances[i],
 DijkstraModes::Undirected);
 }
 
 auto heuristic = [&distances] (node n, unsigned objective) {
 return distances[objective][n];
 };
 
 steady_clock::time_point heuristic_end = steady_clock::now();
 duration<double> heuristic_computation_span
 = duration_cast<duration<double>>(heuristic_end - start);
 
 cout << "Done. (" << heuristic_computation_span.count() << "s)" << endl;
 
 Point absolute_bound(numeric_limits<double>::infinity(), dimension);
 Point relative_bound(numeric_limits<double>::infinity(), dimension);
 
 for(unsigned i = 0; i < dimension; ++i) {
 absolute_bound[i] = numeric_limits<double>::infinity();
 }
 absolute_bound[dimension - 1] = 1.4 * distances[dimension - 1][source];
 
 absolute_bound[0] = 0.01 * absolute_bound[dimension - 1];    // Bündelung
 absolute_bound[1] = 0.01 * absolute_bound[dimension - 1];    // Flughäfen
 absolute_bound[2] = 0.01 * absolute_bound[dimension - 1];    // Tagebau
 absolute_bound[3] = 0.01 * absolute_bound[dimension - 1];    // VSG
 absolute_bound[4] = 0.01 * absolute_bound[dimension - 1];   // Siedlungen
 
 cout << "Running first phase... ";
 
 list<pair<NodeArray<Point *>, NodeArray<edge>>> solutions;
 
 auto callback = [&solutions, &graph, dimension] (NodeArray<Point *>& distances,
 NodeArray<edge>& predecessors) {
 
 NodeArray<Point *> new_distances(graph);
 
 for(auto n : distances.graphOf()->nodes) {
 Point* p = new Point(dimension);
 std::copy(distances[n]->cbegin() + 1, distances[n]->cend(),
 p->begin());
 new_distances[n] = p;
 }
 
 solutions.push_back(make_pair(new_distances, predecessors));
 };
 
 {
 
 EPDualBensonSolver<> weighted_solver;
 
 weighted_solver.Solve(graph, cost_function, source, target, callback);
 
 }
 
 steady_clock::time_point fp_end = steady_clock::now();
 duration<double> fp_computation_span
 = duration_cast<duration<double>>(fp_end - heuristic_end);
 
 cout << "Done. (" << fp_computation_span.count() << "s)" << endl;
 
 vector<list<node>> paths;
 vector<Point> values;
 
 auto path_callback = [&paths] (list<node> path) {
 paths.push_back(path);
 };
 
 auto value_callback = [&values] (Point p) {
 values.push_back(p);
 };
 
 EpSolverMartins solver;
 
 solver.set_path_callback(path_callback);
 solver.set_value_callback(value_callback);
 
 cout << "Running Martins algorithm... ";
 std::flush(cout);
 
 solver.Solve(graph,
 cost_function,
 dimension,
 source,
 target,
 absolute_bound,
 solutions,
 heuristic,
 false);
 
 steady_clock::time_point martins_end = steady_clock::now();
 duration<double> martins_computation_span
 = duration_cast<duration<double>>(martins_end - fp_end);
 
 cout << "Done. (" << martins_computation_span.count() << "s)" << endl;
 
 cout << "Size of the Pareto-frontier: " << solver.solutions().size() << endl;
 
 if(values.size() < 25) {
 for(auto value : values) {
 cout << value << endl;
 }
 }
 
 } else {
 cout << "Unknown algorithm: " << algorithm << endl;
 }
*/
