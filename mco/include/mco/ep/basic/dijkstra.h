/**
 * dijkstra.h
 *
 *  Created on: 25.11.2013
 *      Author: Denis Kurz
 */

#ifndef MCO_DIJKSTRA_H_
#define MCO_DIJKSTRA_H_

#include <functional>
#include <limits>
#include <utility>

#include <ogdf/basic/Graph.h>
#include <ogdf/internal/heap/PairingHeap.h>

#include <mco/basic/point.h>
#include <mco/basic/forward_star.h>

namespace mco {

//! TODO doxygen
template<typename NODE = ogdf::node, typename EDGE = ogdf::edge>
struct DijkstraModes {

    // TODO doxygen
    static std::function<bool(NODE, EDGE)> const Forward;

    // TODO doxygen
    static std::function<bool(NODE, EDGE)> const Backward;

    // TODO doxygen
    static std::function<bool(NODE, EDGE)> const Undirected;

};

template<typename NODE, typename EDGE>
std::function<bool(NODE, EDGE)> const DijkstraModes<NODE, EDGE>::Forward =
[](NODE v, EDGE e) {
    return v == e->source();
};

template<typename NODE, typename EDGE>
std::function<bool(NODE,EDGE)> const DijkstraModes<NODE, EDGE>::Backward =
[](NODE v, EDGE e) {
    return v == e->target();
};
template<typename NODE, typename EDGE>
std::function<bool(NODE,EDGE)> const DijkstraModes<NODE, EDGE>::Undirected =
[](NODE, EDGE) {
    return true;
};

template<typename T, typename C, typename NODE = ogdf::node>
class PairComparator {
public:
    PairComparator& operator=(const PairComparator& that) {
        return *this;
    }

    bool operator()(std::pair<T, NODE>& p1,
                    std::pair<T, NODE>& p2) {
        return comp_(p1.first, p2.first);
    }
private:
    C comp_;
};

//! TODO doxygen
template<typename T, typename C, typename GRAPH = ogdf::Graph, typename NODE = ogdf::node, typename EDGE = ogdf::edge, typename NODE_ARRAY_E = ogdf::NodeArray<ogdf::edge>, typename NODE_ARRAY_T = ogdf::NodeArray<T>>
class Dijkstra {

public:

    //! TODO doxygen
    bool singleSourceShortestPaths(
            GRAPH const &graph,
            std::function<T(EDGE)> weight,
            NODE const source,
            NODE_ARRAY_E &predecessor,
            NODE_ARRAY_T &distance,
            std::function<bool(NODE, EDGE)> mode = DijkstraModes<NODE, EDGE>::Forward) {

        using std::pair;
        using std::make_pair;
        using std::numeric_limits;
        using ogdf::PairingHeap;
        using ogdf::PairingHeapNode;
        using ogdf::NodeArray;

        using queue_element = pair<T, NODE>;

        // Initialization: populate priority queue
        PairingHeap<queue_element, C> queue;
        NodeArray<PairingHeapNode<queue_element>*> qpos(graph);
        // Initialization: set distances
        for(auto v : graph.nodes) {
            distance[v] = numeric_limits<T>::max();
            qpos[v] = queue.push(make_pair(distance[v], v));
        }
        distance[source] = 0;
        queue.decrease(qpos[source], make_pair(distance[source], source));

        // Dijkstra: empty queue, update distances accordingly
        for(int i = 0; i < graph.numberOfNodes(); ++i) {
            auto qe = queue.top();
            queue.pop();
            auto v = qe.second;
            for(auto adj : v->adjEntries) {
                auto e = adj->theEdge();
                if(!mode(v, e)) continue;
                auto w = e->opposite(v);
                T newDist = distance[v] + weight(e);
                if(distance[w] > newDist) {
                    distance[w] = newDist;
                    queue.decrease(qpos[w], make_pair(distance[w], w));
                    predecessor[w] = e;
                }
            }
        }
        
        return true;
    }

};

//! TODO doxygen
class LexDijkstra {

public:

    //! TODO doxygen
    void singleSourceShortestPaths(
            ogdf::Graph const &graph,
            std::function<Point*(ogdf::edge)> weight,
            ogdf::node const source,
            ogdf::NodeArray<Point *> &distance,
            ogdf::NodeArray<ogdf::edge> &predecessor,
            std::function<bool(ogdf::node,ogdf::edge)> mode =DijkstraModes<>::Forward);

};

}

#endif /* MCO_DIJKSTRA_H_ */
