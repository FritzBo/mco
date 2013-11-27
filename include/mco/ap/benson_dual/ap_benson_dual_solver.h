#pragma once
/*
 * ap_benson_dual_solver.h
 *
 *  Created on: 14.10.2013
 *      Author: fritz
 */

#ifndef AP_BENSON_DUAL_SOLVER_H_
#define AP_BENSON_DUAL_SOLVER_H_

#include <memory>
#include <set>
#include <deque>

#include <gurobi_c++.h>

#include <mco/ap/AbstractAPSolver.h>
#include <mco/generic/benson_weightspace/dual_benson_scalarizer.h>

namespace mco {

class APBensonDualSolver : public AbstractAPSolver, public ScalarizationSolver {
public:
	APBensonDualSolver(std::shared_ptr<AssignmentInstance> instance, double epsilon = 1E-6);

	virtual ~APBensonDualSolver();

	void Solve() {
		std::list<Point *> solutions;
		dual_benson_solver_->Calculate_solutions(solutions);

		add_solutions(solutions.begin(), solutions.end());
	}

	double solver_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}

	double vertex_enumeration_time() {
		return dual_benson_solver_->vertex_enumeration_time();
	}

	int number_vertices() {
		return dual_benson_solver_->number_vertices();
	}

	int number_facets() {
		return dual_benson_solver_->number_facets();
	}

protected:
	virtual double Solve_scalarization(Point &weights, Point &value);

private:
	double epsilon_;

	DualBensonScalarizer *dual_benson_solver_;

	ogdf::NodeArray<std::list<ogdf::node> > A_;
	ogdf::NodeArray<ogdf::node> mate_;
	ogdf::NodeArray<ogdf::node> exposed_;
	ogdf::NodeArray<ogdf::node> neighbour_;
	ogdf::NodeArray<ogdf::node> label_;
	ogdf::NodeArray<Point> dual_variables_;
	ogdf::NodeArray<Point> slack_;
	ogdf::NodeArray<unsigned int> count_;

	ogdf::EdgeArray<Point> edge_costs;

	std::deque<ogdf::node> queue_;

	Point null_;
	Point infinity_;

	void augment(ogdf::node);
	bool check_equality_subgraph(Point &cost, Point &dual1, Point &dual2, double epsilon_);
	bool slack_at_least_zero(Point &cost, Point &dual1, Point &dual2, double epsilon_);
	bool better_slack(Point &slack, Point &cost, Point &dual1, Point &dual2, double epsilon_);

	clock_t cycles_;
};

} /* namespace mco */
#endif /* AP_BENSON_DUAL_SOLVER_H_ */
