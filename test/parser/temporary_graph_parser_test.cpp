//
//  lex_hungarian_test.cpp
//  mco
//
//  Created by Fritz Bökler on 04.04.14.
//
//

#include <set>
#include <tuple>
#include <string>

using std::set;
using std::tuple;
using std::string;
using std::get;
using std::make_tuple;

#include <gtest/gtest.h>

using ::testing::Values;

#include <ogdf/basic/Graph.h>

using ogdf::Graph;
using ogdf::node;
using ogdf::EdgeArray;

#include <mco/core/point.h>
#include <mco/benchmarks/temporary_graphs_parser.h>

using mco::Point;
using mco::TemporaryGraphParser;

class TemporaryGraphInstanceFixture
: public ::testing::TestWithParam<tuple<string, unsigned, unsigned>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_nodes_ = get<1>(GetParam());
    const unsigned expected_number_of_edges_ = get<2>(GetParam());
    
};

TEST_P(TemporaryGraphInstanceFixture, NodeEdgeMatch) {
    Graph graph;
    EdgeArray<Point> costs(graph);
    unsigned dimension;
    node source;
    node target;
    
    TemporaryGraphParser parser;
    
    parser.getGraph(filename_,
                    graph,
                    costs,
                    dimension,
                    source,
                    target);
    
    EXPECT_EQ(expected_number_of_nodes_, graph.numberOfNodes());
    EXPECT_EQ(expected_number_of_edges_, graph.numberOfEdges());
}

INSTANTIATE_TEST_CASE_P(TemporaryGraphParsingTest,
                        TemporaryGraphInstanceFixture,
                        Values(
                               make_tuple(string("/Users/fritz/Documents/research/strom/graphs/graph_1000_1000"),
                                          (unsigned) 318,
                                          (unsigned) 1166),
                               make_tuple(string("/Users/fritz/Documents/research/strom/graphs/graph_500_500"),
                                          (unsigned) 1273,
                                          (unsigned) 4880),
                               make_tuple(string("/Users/fritz/Documents/research/strom/graphs/graph_250_250"),
                                          (unsigned) 5121,
                                          (unsigned) 20058)
//                               make_tuple(string("/Users/fritz/Documents/research/strom/graphs/graph_100_100"),
//                                          (unsigned) 83908,
//                                          (unsigned) 202879),
//                               make_tuple(string("/Users/fritz/Documents/research/strom/graphs/graph_50_50"),
//                                          (unsigned) 347360,
//                                          (unsigned) 833615),
//                               make_tuple(string("/Users/fritz/Documents/research/strom/graphs/graph_25_25"),
//                                          (unsigned) 1310225,
//                                          (unsigned) 3235520)
                               ));
