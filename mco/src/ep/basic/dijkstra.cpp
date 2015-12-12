/**
 * dijkstra.cc
 *
 *  Created on: 22.03.2014
 *      Author: Denis Kurz
 */

#include <functional>
#include <limits>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/basic/lex_point_comparator.h>
#include <mco/ep/basic/dijkstra.h>

using std::function;
using std::numeric_limits;
using std::make_pair;

using ogdf::Graph;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using ogdf::PairingHeap;
using ogdf::PairingHeapNode;

namespace mco {

void LexDijkstra::singleSourceShortestPaths(
        Graph const &graph,
        function<Point*(edge)> weight,
        node const source,
        NodeArray<Point *> &distance,
        NodeArray<edge> &predecessor,
        function<bool(node,edge)> mode) {

    using queue_element = std::pair<Point *, node>;
    
    LexPointComparator less;
    
    unsigned dim = weight(graph.chooseEdge())->dimension();

    PairingHeap<queue_element, PairComparator<Point*, LexPointComparator>> queue;
    NodeArray<PairingHeapNode<queue_element>*> qpos(graph);
    for(auto v : graph.nodes) {
        for(unsigned i = 0; i < dim; ++i) {
            (*distance[v])[i] = numeric_limits<double>::max();
        }
        predecessor[v] = nullptr;
        qpos[v] = queue.push(make_pair(distance[v], v));
    }

    *distance[source] = Point(0.0, distance[source]->dimension());
    
    queue.decrease(qpos[source], make_pair(distance[source], source));

    Point tmp(dim);
    for(int i = 0; i < graph.numberOfNodes(); ++i) {
        queue_element qe = queue.top();
        auto v = qe.second;
        queue.pop();
        for(auto adj : v->adjEntries) {
            edge e = adj->theEdge();
            if(!mode(v, e)) continue;
            node w = e->opposite(v);
            tmp = *distance[v];
            tmp += *weight(e);
            if(less(&tmp, distance[w])) {
                *distance[w] = tmp;
                queue.decrease(qpos[w], make_pair(distance[w], w));
                predecessor[w] = e;
            }
        }
    }
}

function<bool(node,edge)> const DijkstraModes::Forward =
    [](node v, edge e) {
        return v == e->source();
    };

function<bool(node,edge)> const DijkstraModes::Backward =
    [](node v, edge e) {
        return v == e->target();
    };

function<bool(node,edge)> const DijkstraModes::Undirected =
    [](node, edge) {
        return true;
    };

}

