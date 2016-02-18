//
//  forward_star.cpp
//  mco
//
//  Created by Fritz Bökler on 09.01.16.
//
//

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

using std::ifstream;
using std::string;
using std::map;
using std::vector;
using std::pair;
using std::make_pair;

#include <mco/basic/forward_star.h>

namespace mco
{

ForwardStar::ForwardStar(ogdf::Graph& ogdf_graph, bool directed)
:   no_nodes(0),
    no_edges(0),
    nodes(*this),
    heads_(*this),
    tails_(*this)
{
    no_edges = ogdf_graph.numberOfEdges();
    no_nodes = ogdf_graph.numberOfNodes();

    first_edge_.reserve(no_nodes);
    heads_.objects_.resize(no_edges);
    tails_.objects_.resize(no_edges);

    int max_node_id = 0;
    for(auto n : ogdf_graph.nodes)
    {
        max_node_id = std::max(max_node_id, n->index());
    }

    vector<node> node_ids(max_node_id);

    unsigned i = 0;
    for(auto n : ogdf_graph.nodes)
    {
        node_ids[n->index()] = i;
        ++i;
    }

    unsigned j = 0;
    for(auto n: ogdf_graph.nodes)
    {
        first_edge_[node_ids[n->index()]] = j;

        for(auto adj : n->adjEntries)
        {
            auto e = adj->theEdge();
            if(e->target() == n && directed)
            {
                continue;
            }

            heads_.objects_[j] = node_ids[e->opposite(n)->index()];
            tails_.objects_[j] = node_ids[n->index()];

            ++j;
        }
    }

    first_edge_[no_nodes] = no_edges;
}

void ForwardStarFileReader::read(std::string filename,
                                 ForwardStar& graph,
                                 FSNodeArray<int>& extern_node_ids,
                                 FSEdgeArray<Point>& weights,
                                 unsigned& dimension,
                                 node& source,
                                 node& target,
                                 bool directed)
{

    ifstream file(filename);

    if(!file.good()) {
#ifndef NDEBUG
        std::cerr << "Could not open file " << filename << std::endl;
#endif
        throw string("Could not open file ") + filename;
    }

    unsigned num_nodes, num_edges, source_id, target_id;

    file >> num_nodes; //number of nodes
    file >> num_edges; //number of edges
    file >> dimension;
    file >> source_id;
    file >> target_id;

    int node1_ref, node2_ref;

    map<int, node> node_ids;
    vector<vector<pair<node, Point>>> neighbors(num_nodes);
    extern_node_ids.objects_.resize(num_nodes);

    node id = 0;
    for(unsigned i = 0; i < num_edges; ++i ) {
        file >> node1_ref >> node2_ref;

        if(node_ids.count(node1_ref) == 0)
        {
            node_ids[node1_ref] = id;
            extern_node_ids[node_ids[node1_ref]] = node1_ref;
            ++id;
        }

        if(node_ids.count(node2_ref) == 0)
        {
            node_ids[node2_ref] = id;
            extern_node_ids[node_ids[node2_ref]] = node2_ref;
            ++id;
        }

        // FIXME
        Point new_point(dimension);
        for(unsigned j = 0; j < dimension; ++j) {
            file >> new_point[j];
        }

        if(!directed)
        {
            neighbors[node_ids[node2_ref]].push_back(make_pair(node_ids[node1_ref], new_point));
        }
        neighbors[node_ids[node1_ref]].push_back(make_pair(node_ids[node2_ref], std::move(new_point)));
    }

    if(!directed)
    {
        num_edges *= 2;
    }

    source = node_ids[source_id];
    target = node_ids[target_id];

    //close the file
    file.close();

    // Create forward star representation
    graph.no_nodes = num_nodes;
    graph.no_edges = num_edges;

    graph.first_edge_.resize(num_nodes + 1);
    graph.heads_.objects_.resize(num_edges);
    graph.tails_.objects_.resize(num_edges);
    weights.objects_.resize(num_edges);

    unsigned j = 0;
    for(unsigned i = 0; i < num_nodes; ++i)
    {
        graph.first_edge_[i] = j;

        for(unsigned k = 0; k < neighbors[i].size(); ++k)
        {
            assert(j < num_edges);

            graph.heads_[j] = neighbors[i][k].first;
            graph.tails_[j] = i;
            weights[j] = std::move(neighbors[i][k].second);

            ++j;
        }
    }

    graph.first_edge_[num_nodes] = num_edges;

}

}