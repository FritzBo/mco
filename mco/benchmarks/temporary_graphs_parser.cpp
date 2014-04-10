//
//  temporary_graphs_parser.cpp
//  mco
//
//  Created by Fritz Bökler on 08.04.14.
//
//

#include <mco/benchmarks/temporary_graphs_parser.h>

#include <map>

using std::map;
using std::pair;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::EdgeArray;

#include <mco/basic/point.h>

using mco::Point;

namespace mco {
    
void TemporaryGraphParser::getGraph(string filename,
                                    Graph &graph,
                                    EdgeArray<Point> &weights,
                                    unsigned& dimension,
                                    ogdf::node& source,
                                    ogdf::node& target) {
    
    ifstream file(filename);
    
    if(!file.good()) {
        cerr << "Could not open file " << filename_ << endl;
        throw string("Could not open file ") + filename_;
    }
    
    unsigned num_nodes, num_edges, source_id, target_id;
    
    file >> num_nodes; //number of nodes
    file >> num_edges; //number of variables
    file >> dimension;
    file >> source_id;
    file >> target_id;
    
    int node1_ref, node2_ref;
    node node1, node2;
    edge e;
    map<int, node> nodes_added;
    
    for( int i = 0; i < num_edges; ++i ) {
        file >> node1_ref >> node2_ref;
        
        if(nodes_added.count(node1_ref) == 0) {
            node1 = graph.newNode(node1_ref);
            nodes_added.insert(pair<int, node>(node1_ref, node1));
        } else
            node1 = nodes_added[node1_ref];
        
        if(nodes_added.count(node2_ref) == 0) {
            node2 = graph.newNode(node2_ref);
            nodes_added.insert(pair<int, node>(node2_ref, node2));
        } else
            node2 = nodes_added[node2_ref];
        
        // FIXME
        Point new_point(dimension + 1);
        file >> new_point[0];
        new_point[0] += 1;
        
        for(int j = 1; j < dimension; j++) {
            file >> new_point[j];
        }
        
        new_point[dimension] = 1;
        
        e = graph.newEdge(node1, node2);
        weights[e] = std::move(new_point);
    }
    
    source = nodes_added[source_id];
    target = nodes_added[target_id];
    
    //FIXME
    ++dimension;
    
    //close the file
    file.close();
}
    
}