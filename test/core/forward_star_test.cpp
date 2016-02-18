//
//  forward_star_test.cpp
//  mco
//
//  Created by Fritz Bökler on 10.01.16.
//
//

#include <limits>
#include <utility>

using std::numeric_limits;
using std::make_pair;

#include <gtest/gtest.h>

#include <ogdf/internal/heap/PairingHeap.h>

using ogdf::PairingHeap;
using ogdf::PairingHeapNode;

#include <mco/basic/point.h>
#include <mco/basic/forward_star.h>

using mco::Point;
using mco::LexPointComparator;
using mco::ForwardStar;
using mco::node;
using mco::edge;
using mco::FSEdgeArray;
using mco::FSNodeArray;
using mco::ForwardStarFileReader;

TEST(ForwardStarTest, FileReadTest)
{
    ForwardStarFileReader reader;

    ForwardStar graph;
    FSNodeArray<int> extern_ids(graph);
    FSEdgeArray<Point> weights(graph);
    unsigned dimension;
    node source, target;

    reader.read("../instances/ep/graph_1000_1000",
                graph,
                extern_ids,
                weights,
                dimension,
                source,
                target,
                false);

    EXPECT_EQ(graph.numberOfNodes(), 318);
    EXPECT_EQ(graph.numberOfEdges(), 1166);
    EXPECT_EQ(source, 190);
    EXPECT_EQ(target, 35);

    unsigned i = 0;
    for(edge e : graph.adj_edges(190))
    {
        EXPECT_TRUE(e == 695 ||
                    e == 696 ||
                    e == 697);
        ++i;
    }
    EXPECT_EQ(i, 3);

    i = 0;
    for(edge e : graph.adj_edges(100))
    {
        EXPECT_TRUE(e == 367 ||
                    e == 368 ||
                    e == 369 ||
                    e == 370);

        ++i;
    }
    EXPECT_EQ(i, 4);

}

template<typename T, typename C>
class PairComparator {
public:
    PairComparator& operator=(const PairComparator& that) {
        return *this;
    }

    bool operator()(std::pair<T, node>& p1,
                    std::pair<T, node>& p2) {
        return comp_(p1.first, p2.first);
//        return comp_(p2.first, p1.first);
    }
private:
    C comp_;
};


TEST(ForwardStarTest, DijkstraTest)
{
    ForwardStarFileReader reader;

    ForwardStar graph;
    FSNodeArray<int> extern_ids(graph);
    FSEdgeArray<Point> weights(graph);
    unsigned dimension;
    node source, target;

    reader.read("../instances/ep/grid80_100_4",
                graph,
                extern_ids,
                weights,
                dimension,
                source,
                target,
                true);

    FSNodeArray<node> predecessor(graph);
    FSNodeArray<Point*> distance(graph);

    for(auto n : graph.nodes)
    {
        distance[n] = new Point(dimension);
    }

    using queue_element = std::pair<Point *, node>;

    LexPointComparator less;

    PairingHeap<queue_element, PairComparator<Point*, LexPointComparator>> queue;
    FSNodeArray<PairingHeapNode<queue_element>*> qpos(graph);

    for(auto v : graph.nodes) {
        for(unsigned i = 0; i < dimension; ++i) {
            (*distance[v])[i] = numeric_limits<double>::max();
        }
        predecessor[v] = source;
        qpos[v] = queue.push(make_pair(distance[v], v));
    }

    *distance[source] = Point(0.0, distance[source]->dimension());

    queue.decrease(qpos[source], make_pair(distance[source], source));

    Point tmp(dimension);
    for(unsigned i = 0; i < graph.numberOfNodes(); ++i) {
        queue_element qe = queue.top();
        auto v = qe.second;
        queue.pop();
        for(edge e : graph.adj_edges(v)) {
            if(graph.head(e) == v)
            {
                continue;
            }
            node w = graph.head(e);
            tmp = *distance[v];
            tmp += weights[e];
            if(less(&tmp, distance[w])) {
                *distance[w] = tmp;
                queue.decrease(qpos[w], make_pair(distance[w], w));
                predecessor[w] = e;
            }
        }
    }

    EXPECT_EQ(round((*distance[target])[0]), 2656);
    EXPECT_EQ(round((*distance[target])[1]), 5313);

//    std::cout << *distance[target] << std::endl;
}