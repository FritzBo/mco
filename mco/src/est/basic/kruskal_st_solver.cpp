/*
 * prim_st_solver.cpp
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#include <mco/est/basic/kruskal_st_solver.h>

#include <vector>
#include <list>

using std::vector;
using std::list;

using ogdf::edge;
using ogdf::node;

namespace mco {

LexKruskalSTSolver::~LexKruskalSTSolver() {
}

void LexKruskalSTSolver::
Solve(const ogdf::Graph & g,
      std::function<Point*(ogdf::edge)> costs,
      unsigned dimension) {

	vector<edge> sorted_edges;

	cost_ = Point(dimension);
	spanning_tree_.clear();
	spanning_tree_.resize(g.numberOfNodes() - 1);

    sets_.init(g);

	for(auto e : g.edges) {
		sorted_edges.push_back(e);
	}

    auto comparator = [&costs] (edge e1, edge e2) {
        return LexPointComparator::is_lex_le(costs(e1),
                                             costs(e2),
                                             0.0);
    };

    stable_sort(sorted_edges.begin(),
                sorted_edges.end(),
                comparator);

	for(auto n : g.nodes) {
		make_set(n);
	}

    int i = 0;
	for(edge e: sorted_edges) {

        if(i > g.numberOfNodes()) {
            break;
        }

		if(find_set(e->source()) != find_set(e->target())) {
			spanning_tree_.push_back(e);
			set_union(e->source(), e->target());
			cost_ += (*costs(e));
            ++i;
		}
	}
}

const vector<edge> & LexKruskalSTSolver::spanning_tree() {
	return spanning_tree_;
}

Point LexKruskalSTSolver::min_cost() {
	return cost_;
}

void LexKruskalSTSolver::make_set(ogdf::node n) {
    sets_[n] = n;
}

ogdf::node LexKruskalSTSolver::find_set(ogdf::node n) {
	node parent = n;
	list<node> path_nodes;

    // Path compression
	while(sets_[parent] != parent) {
		path_nodes.push_back(parent);
		parent = sets_[parent];
    }

	for(node n: path_nodes) {
		sets_[n] = parent;
    }

	return parent;
}

void LexKruskalSTSolver::set_union(ogdf::node u, ogdf::node v) {
	sets_[find_set(v)] = sets_[find_set(u)];
}

} /* namespace mco */
