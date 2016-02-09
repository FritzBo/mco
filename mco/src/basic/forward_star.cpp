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

void ForwardStarFileReader::read(std::string filename,
                                 ForwardStar& graph,
                                 FSEdgeArray<Point>& weights,
                                 unsigned& dimension,
                                 node& source,
                                 node& target)
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
    file >> num_edges; //number of variables
    file >> dimension;
    file >> source_id;
    file >> target_id;

    int node1_ref, node2_ref;
    vector<vector<pair<node, Point>>> neighbors(num_nodes);

    for(unsigned i = 0; i < num_edges; ++i ) {
        file >> node1_ref >> node2_ref;

        // FIXME
        Point new_point(dimension);
        for(unsigned j = 0; j < dimension; ++j) {
            file >> new_point[j];
        }

        neighbors[node1_ref].push_back(make_pair(node2_ref, std::move(new_point)));
    }

    source = source_id;
    target = target_id;

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
            graph.heads_[j] = neighbors[i][k].first;
            graph.tails_[j] = i;
            weights[j] = std::move(neighbors[i][k].second);

            ++j;
        }
    }

    graph.first_edge_[num_nodes] = num_edges;

}

}