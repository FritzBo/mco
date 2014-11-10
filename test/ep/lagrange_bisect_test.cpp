//
//  ep_benson_dual_test.cpp
//  mco
//
//  Created by Fritz Bökler on 08.08.2014
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
using ogdf::edge;
using ogdf::EdgeArray;

#include <mco/basic/point.h>
#include <mco/benchmarks/temporary_graphs_parser.h>

#include <mco/ep/basic/ep_lagrange_bisect.h>

using mco::Point;
using mco::TemporaryGraphParser;

using mco::EpLagrangeBisect;

class BoundedParetoInstanceTestFixture
: public ::testing::TestWithParam<tuple<string, unsigned, Point>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_minimizers_ = get<1>(GetParam());
    const Point bounds_ = get<2>(GetParam());
};

TEST_P(BoundedParetoInstanceTestFixture, CountMatch) {
    Graph graph;
    EdgeArray<Point> costs(graph);
    unsigned dimension;
    node source;
    node target;
    
    TemporaryGraphParser parser;
    
    parser.getGraph(filename_, graph, costs, dimension, source, target);
    
    EpLagrangeBisect bisect;
    
    auto weight_function = [costs] (edge e) -> const Point& {
        return costs(e);
    };
    
    Point lambda(1/(double) dimension, dimension);
    
    double value = bisect.find_lagrange_multi(graph,
                                              weight_function,
                                              dimension,
                                              source,
                                              target,
                                              true,
                                              bounds_,
                                              lambda);
    
    
#ifndef NDEBUG
    cout << "Approximate Lagrange Multiplier: " << lambda << endl;
    cout << "Value: " << value << endl;
#endif
    
}

INSTANTIATE_TEST_CASE_P(InstanceTests,
                        BoundedParetoInstanceTestFixture,
                        Values(
                               //                               make_tuple(string("../../../instances/ep/graph_1000_1000"), (unsigned) 1421),
                               make_tuple(string("../../../instances/ep/grid50_50_7"), (unsigned) 2,
                                          Point({251900, 237542}))
//                               make_tuple(string("../../../instances/ep/grid50_50_7"), (unsigned) 13)
                               ));
