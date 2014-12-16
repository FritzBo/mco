#pragma once
/*
 * prim_st_solver.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef KRUSKAL_ST_SOLVER_H_
#define KRUSKAL_ST_SOLVER_H_

#include <map>
#include <vector>
#include <functional>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class LexKruskalSTSolver {
public:
	virtual ~LexKruskalSTSolver();

	void Solve(const ogdf::Graph & g,
               std::function<Point*(ogdf::edge)> costs,
               unsigned dimension);

	const std::vector<ogdf::edge> & spanning_tree();
	Point min_cost();

private:
	std::vector<ogdf::edge> spanning_tree_;
	Point cost_;

	ogdf::NodeArray<ogdf::node> sets_;
	void make_set(ogdf::node n);
	ogdf::node find_set(ogdf::node n);
	void set_union(ogdf::node u, ogdf::node v);
};

} /* namespace mco */
#endif /* PRIM_ST_SOLVER_H_ */
