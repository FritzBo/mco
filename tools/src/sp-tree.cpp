//
//  main.cpp
//  cli
//
//  Created by Fritz Bökler on 25.04.2016.
//
//

#include <thread>
#include <string>
#include <list>
#include <vector>

using std::string;
using std::list;
using std::vector;

#include <tclap/CmdLine.h>

using TCLAP::CmdLine;
using TCLAP::UnlabeledValueArg;
using TCLAP::SwitchArg;
using TCLAP::ValueArg;
using TCLAP::ArgException;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::EdgeArray;
using ogdf::NodeArray;
using ogdf::node;
using ogdf::edge;

#include <mco/basic/point.h>
#include <mco/ep/basic/dijkstra.h>
#include <mco/benchmarks/temporary_graphs_parser.h>
#include <mco/cli/parse_util.h>

using mco::Point;
using mco::LexDijkstra;
using mco::DijkstraModes;
using mco::TemporaryGraphParser;
using mco::LexPointComparator;


using cost_function_type = std::function<Point*(edge)>;

void print_sidetracks(std::vector<edge>& sidetracks,
                      cost_function_type costs,
                      bool weighted)
{
    for(auto e : sidetracks)
    {
        unsigned start = weighted ? 1 : 0;
        assert(costs(e)->dimension() > 0);
        for(unsigned i = start; i < costs(e)->dimension(); ++i)
        {
            std::cout << (*costs(e))[i] << ", ";
        }
        std::cout << e->source() << ", " << e->target() << ", " << std::endl;
    }
}

void find_sidetracks(Graph& graph,
                     cost_function_type costs,
                     NodeArray<edge>& predecessors,
                     NodeArray<Point*>& distances,
                     vector<edge>& sidetracks)
{
    EdgeArray<bool> sp_edges(graph, false);

    for(auto n : graph.nodes)
    {
        if(predecessors[n] != nullptr)
        {
            sp_edges[predecessors[n]] = true;
        }
    }

    Point epsilon(1E-6, costs(graph.chooseEdge())->dimension());
    LexPointComparator comp;

    for(auto e : graph.edges)
    {
        if(sp_edges[e])
        {
            continue;
        }

        node tail = e->source();
        node head = e->target();

        if(comp(*distances[tail] - *distances[head] + *costs(e), epsilon))
        {
            sidetracks.push_back(e);
        }
    }
}

void print_root_to_leaf_path(NodeArray<edge>& predecessors,
                             NodeArray<Point*>& distances,
                             node source,
                             node n,
                             bool weighted)
{
    unsigned start = weighted ? 1 : 0;
    assert(distances[n]->dimension() > 0);
    for(unsigned i = start; i < distances[n]->dimension(); ++i) {
        std::cout << (*distances[n])[i] << ", ";
    }

    list<node> path;
    auto curr = n;
    path.push_front(curr);

    while(curr != source) {

        curr = predecessors[curr]->opposite(curr);
        path.push_front(curr);
    }

    for(auto n : path) {
        std::cout << n << ", ";
    }

    std::cout << std::endl;
}

void compute_leaves(Graph& graph,
                    NodeArray<edge>& predecessors,
                    vector<node>& leaves)
{
    NodeArray<bool> leaf_candidate(graph, true);
    for(auto n : graph.nodes)
    {
        if(predecessors[n] != nullptr)
        {
            leaf_candidate[predecessors[n]->opposite(n)] = false;
        }
    }
    for(auto n: graph.nodes)
    {
        if(leaf_candidate[n])
        {
            leaves.push_back(n);
        }
    }
}

Point parse_weight(ValueArg<string>& argument,
                   unsigned dimension) {

    vector<string> tokens;
    mco::tokenize(argument.getValue(), tokens, ":");

    if(tokens.size() != dimension) {
        std::cout << "Error: Wrong number of weights" << std::endl;
        return Point();
    }

    Point r(dimension);

    for(unsigned i = 0; i < dimension; ++i) {
        r[i] = stod(tokens[i]);
    }
    
    return r;
}

int main(int argc, char** argv)
{
    // todo: remove
#ifndef NDEBUG
    std::this_thread::sleep_for(std::chrono::seconds(1));
#endif

    try {
        CmdLine cmd("Computes a lexicographic minimal shortest path tree with optional weights.", ' ', "0.1");

        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        SwitchArg is_directed_arg("d", "directed", "Should the input be interpreted as a directed graph?", false);

        SwitchArg count_only_arg("c", "count", "Should only the number of leaves be counted?", false);

        ValueArg<string> weight_arg("W", "weight", "Adding a weighted-sum as first objective", false, "", "weight");


        cmd.add(file_name_argument);
        cmd.add(is_directed_arg);
        cmd.add(weight_arg);
        cmd.add(count_only_arg);

        cmd.parse(argc, argv);

        string file_name = file_name_argument.getValue();
        bool is_directed = is_directed_arg.getValue();
        bool count_only = count_only_arg.getValue();

        Graph graph;
        EdgeArray<Point> costs(graph);
        unsigned dimension;
        node source, target;

        TemporaryGraphParser parser;

        parser.getGraph(file_name, graph, costs, dimension, source, target);

        NodeArray<Point*> distances(graph);
        NodeArray<edge> predecessors(graph);

        LexDijkstra sssp;

        std::function<Point*(edge)> cost_function;
        EdgeArray<Point>* weighted_cost = nullptr;
        bool weighted;

        if(weight_arg.isSet()) {

            Point weight = parse_weight(weight_arg, dimension);

            weighted_cost = new ogdf::EdgeArray<Point>(graph);

            for(auto e : graph.edges) {
                (*weighted_cost)[e] = Point(dimension + 1);

                (*weighted_cost)[e][0] = weight * costs[e];

                for(unsigned i = 1; i <= dimension; ++i) {
                    (*weighted_cost)[e][i] = costs[e][i - 1];
                }
            }

            for(auto n : graph.nodes) {
                distances[n] = new Point(dimension + 1);
            }

            cost_function = [weighted_cost] (edge e) {
                return &(*weighted_cost)[e];
            };

            sssp.singleSourceShortestPaths(graph,
                                           cost_function,
                                           source,
                                           distances,
                                           predecessors,
                                           is_directed ? DijkstraModes::Forward : DijkstraModes::Undirected);

            weighted = true;



        } else {
            cost_function = [&costs] (edge e) {
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

            weighted = false;
            
        }

        vector<node> leaves;
        compute_leaves(graph, predecessors, leaves);

        if(count_only)
        {
            std::cout << leaves.size() << std::endl;
        }
        else
        {

            print_root_to_leaf_path(predecessors, distances, source, target, weighted);

            for(auto n : leaves)
            {
                if(n != source)
                {
                    print_root_to_leaf_path(predecessors, distances, source, n, weighted);
                }
            }

            vector<edge> sidetracks;
            find_sidetracks(graph, cost_function, predecessors, distances, sidetracks);

            print_sidetracks(sidetracks, cost_function, weighted);
        }


        for(auto n : graph.nodes) {
            delete distances[n];
        }

        delete weighted_cost;
        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    return 0;
}