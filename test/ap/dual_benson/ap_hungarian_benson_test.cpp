//
//  ap_hungarian_benson_test.cpp
//  mco
//
//  Created by Fritz Bökler on 06.04.14.
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

#include <mco/basic/point.h>
#include <mco/ap/basic/ap_instance.h>
#include <mco/ap/benson_dual/ap_benson_dual_solver.h>
#include <mco/generic/benson_dual/ove_cdd.h>
#include <mco/benchmarks/mcap_parser.h>

#include <mco/ap/brute_force/ap_brute_force_solver.h>

using mco::Point;
using mco::MCAPParser;
using mco::AssignmentInstance;
using mco::APBensonDualSolver;
using mco::OnlineVertexEnumeratorCDD;
using mco::EqualityPointComparator;

using mco::APBruteForceSolver;

using mco::EqualityPointComparator;

class ParetoInstanceTestFixture
: public ::testing::TestWithParam<tuple<string, unsigned>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_minimizers_ = get<1>(GetParam());
    
};

TEST_P(ParetoInstanceTestFixture, CountMatch) {
    Graph graph;
    EdgeArray<Point*> costs(graph);
    set<node> agents;
    
    MCAPParser parser(filename_);
    
    AssignmentInstance instance = parser.get_instance(graph,
                                                      costs,
                                                      agents);

    APBruteForceSolver solver2(instance);
    solver2.Solve();

//    for(auto s : solver2.solutions()) {
//        cout << "1\t" << s.second[0];
//        for(unsigned i = 1; i < instance.dimension(); ++i) {
//            cout << "\t" << s.second[i];
//        }
//        cout << endl;
//    }
//
//    cout << solver2.solutions().size() << endl;

    APBensonDualSolver<OnlineVertexEnumeratorCDD> solver;
    
    solver.Solve(instance);

    cout << solver.number_facets() << endl;
    
    EXPECT_EQ(expected_number_of_minimizers_, solver.solutions().size());

    cout << endl;

    for(auto s : solver.solutions()) {
        bool found = false;
        for(auto s2 : solver2.solutions()) {
            if(EqualityPointComparator(1E-5)(s.second, s2.second)) {
                found = true;
                break;
            }
        }
        if(found) {
            break;
        } else {
            cout << "Not found: " << s.second << endl;
            ASSERT_TRUE(false);
        }
        // (26, 90, 75, 100, 99, 78)
    }

    for(auto e : graph.edges) {
        delete costs(e);
    }
}

INSTANTIATE_TEST_CASE_P(InstanceTests,
                        ParetoInstanceTestFixture,
                        Values(
                               make_tuple(string("../../../instances/ap/ap_6_8_6"), 295),
                               make_tuple(string("../../../instances/ap/pge_3_10_2"), 28),
                               make_tuple(string("../../../instances/ap/pge_3_05_5"), 12),
                               make_tuple(string("../../../instances/ap/map_1_2_0"), 1)
                               ));
