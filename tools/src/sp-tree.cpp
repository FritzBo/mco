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
#include <deque>

using std::string;
using std::list;
using std::vector;
using std::deque;

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
                      EdgeArray<bool>& filter,
                      cost_function_type costs,
                      bool weighted)
{
    for(auto e : sidetracks)
    {
        if(!filter[e])
        {
            continue;
        }

        unsigned start = weighted ? 1 : 0;
        assert(costs(e)->dimension() > 0);
        for(unsigned i = start; i < costs(e)->dimension(); ++i)
        {
            std::cout << (*costs(e))[i] << ", ";
        }
        std::cout << e->source() << ", " << e->target() << ", " << std::endl;
    }
}

void filter_edges(Graph& graph,
                  node target,
                  NodeArray<edge>& predecessors,
                  vector<edge>& sidetracks,
                  EdgeArray<bool>& filter)
{

    // Construct directed graph

    NodeArray<vector<edge>> back_edges(graph);

    for(auto n : graph.nodes)
    {
        if(predecessors[n] != nullptr)
        {
            back_edges[n].push_back(predecessors[n]);
        }
    }

    for(auto e : sidetracks)
    {
        back_edges[e->target()].push_back(e);
    }

    // BFS
    deque<node> queue;
    NodeArray<bool> seen(graph, false);
    queue.push_back(target);
    node curr = nullptr;

    while(!queue.empty())
    {
        curr = queue.front();
        queue.pop_front();
        if(seen[curr])
        {
            continue;
        }
        seen[curr] = true;

        for(auto e : back_edges[curr])
        {
            if(e->source() == curr)
            {
                continue;
            }

            filter[e] = true;

            if(!seen[e->source()])
            {
                queue.push_back(e->source());
            }
        }
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

    mco::EqualityPointComparator eq(1E-6);

    for(auto e : graph.edges)
    {
        if(sp_edges[e])
        {
            continue;
        }

        node tail = e->source();
        node head = e->target();

        if(eq(*distances[tail], *distances[head] + *costs(e)))
        {
            sidetracks.push_back(e);
        }
    }
}

void print_root_to_leaf_path(EdgeArray <Point>& costs,
                             NodeArray<edge>& predecessors,
                             NodeArray<Point*>& distances,
                             NodeArray<bool>& printed_nodes,
                             EdgeArray<bool> filter,
                             node source,
                             node n,
                             unsigned dimension,
                             bool weighted)
{
    Point cost(dimension);
    deque<node> path;
    auto curr = n;
    path.push_front(curr);

    while(curr != source&& !printed_nodes[curr]) {

        cost += costs[predecessors[curr]];

        if(!filter[predecessors[curr]])
        {
            path.clear();

            for(unsigned i = 0; i < cost.dimension(); ++i)
            {
                cost[i] = 0;
            }
        }

        curr = predecessors[curr]->opposite(curr);
        path.push_front(curr);
    }

    if(path.size() > 1)
    {
        for(unsigned i = 0; i < cost.dimension(); ++i) {
            std::cout << cost[i] << ", ";
        }

        for(auto n : path) {
            std::cout << n << ", ";
            printed_nodes[n] = true;
        }

        std::cout << std::endl;
    }
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

void double_edges(Graph& graph,
                  EdgeArray<Point>& costs)
{
    vector<edge> edges;
    for(auto e : graph.edges)
    {
        edges.push_back(e);
    }
    for(auto e : edges)
    {
        auto new_edge = graph.newEdge(e->target(), e->source());
        costs[new_edge] = costs[e];
    }
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

        if(!is_directed)
        {
            double_edges(graph, costs);
        }

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
                                           DijkstraModes::Forward);

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
                                           DijkstraModes::Forward);

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
            vector<edge> sidetracks;
            find_sidetracks(graph, cost_function, predecessors, distances, sidetracks);

            EdgeArray<bool> filter(graph, false);
            filter_edges(graph, target, predecessors, sidetracks, filter);

            NodeArray<bool> printed_nodes(graph, false);
            print_root_to_leaf_path(costs, predecessors, distances, printed_nodes, filter, source, target, dimension, weighted);

            for(auto n : leaves)
            {
                if(n != source)
                {
                    print_root_to_leaf_path(costs, predecessors, distances, printed_nodes, filter, source, n, dimension, weighted);
                }
            }

            print_sidetracks(sidetracks, filter, cost_function, weighted);
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