/*
 * ap_brute_force_solver.cpp
 *
 *  Created on: 06.11.2013
 *      Author: fritz
 */

#include <mco/ap/brute_force/ap_brute_force_solver.h>

#include <set>
#include <list>

using std::set;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::node;

namespace mco {

void APBruteForceSolver::Solve() {
	agent_list.clear();
	for(auto agent: *instance().agents())
		agent_list.push_back(agent);

	set<node> jobs;
	list<edge> current_matching;
	Point current_cost(instance().dimension());
	list<list<edge>> matchings;
	list<Point> costs;

	recursive_find(0, jobs, matchings, costs, current_matching, current_cost);

	for(Point cost : costs)
		add_solution(new Point(cost));
}

void APBruteForceSolver::recursive_find(unsigned int agent_index, set<node> &jobs, list<list<edge>> &matchings, list<Point> &costs, list<edge> &current_matching, Point &current_cost) {
//	cout << agent_index << ", " << jobs.size() << ", " << current_cost << ", " << matchings.size() << endl;
	edge e;
	forall_adj_edges(e, agent_list[agent_index]) {
		if(jobs.count(e->target()) > 0)
			continue;

		jobs.insert(e->target());
		current_matching.push_back(e);
		current_cost += *instance().weights()->operator [](e);

		if(agent_index == instance().agents()->size() - 1) {
			bool dominated = false;

			auto cost_it = costs.begin();
			auto matching_it = matchings.begin();
			while(cost_it != costs.end()) {
				if(cost_it->Dominates(current_cost, 1E-9)) {
					dominated = true;
					break;
				}
				if(current_cost.Dominates(*cost_it, 1E-9)) {
					costs.erase(cost_it);
//					matchings.erase(matching_it);
				}

				cost_it++;
				matching_it++;
			}

			if(dominated) {
				current_matching.pop_back();
				current_cost -= *instance().weights()->operator [](e);
				jobs.erase(e->target());
				break;
			}

			matchings.push_back(current_matching);
			costs.push_back(current_cost);
		} else {

			recursive_find(agent_index + 1,
					jobs, matchings, costs,
					current_matching, current_cost);
		}

		current_matching.pop_back();
		current_cost -= *instance().weights()->operator [](e);
		jobs.erase(e->target());
	}
}
} /* namespace mco */
