/*
 * martins.cpp
 *
 *  Created on: 23.03.2013
 *      Author: fritz
 */

#include <mco/ep/martins/martins.h>

#include <queue>
#include <vector>
#include <set>
#include <list>

using std::priority_queue;
using std::vector;
using std::set;
using std::list;

#include <ogdf/basic/Graph.h>

using ogdf::AdjElement;
using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;

#include <mco/ep/basic/ep_instance.h>
#include <mco/core/point.h>
#include <mco/ep/martins/label.h>

namespace mco {

// TODO: Functor or lambda
bool lexicographic_smaller_label(const Label *p1, const Label *p2) {
    // return *p1->point < *p2->point;
    return LexPointComparator::is_lex_le(*p1->point, *p2->point, 0);
}

using LabelPriorityQueue = priority_queue<const Label *, vector<const Label *>, decltype(&lexicographic_smaller_label)>;

// TODO: FIXME
void EpSolverMartins::Solve() {
	LabelPriorityQueue lex_min_label(&lexicographic_smaller_label);
	NodeArray<list<Label *>> labels(instance().graph());

	set<const Label *> labels_in_queue;
    
    ComponentwisePointComparator comp_leq(epsilon_, false);

	Label *null_label = new Label(Point::Null(instance().dimension()), instance().source(), nullptr);
	labels[instance().source()].push_back(null_label);
	lex_min_label.push(null_label);
	labels_in_queue.insert(null_label);

	while(!lex_min_label.empty()) {
		const Label *label = lex_min_label.top();
		lex_min_label.pop();

		if(label->mark_dominated) {
			delete label;
			continue;
		}

		const Point *label_cost = label->point;
		node n = label->n;

//		cout << endl << n << ": ";

		AdjElement *adj;
		forall_adj(adj, n) {
			edge e = adj->theEdge();

			if(e->isSelfLoop())
				continue;

			node v = e->target();

			if(v == n)
				continue;

			if(v == instance().source())
				continue;

//			cout << v << ", ";

			const Point *edge_cost = instance().weights()(e);
			const Point *new_cost = new Point(*edge_cost + *label_cost);			// Owner

			list<Label *> nondominatedLabels;
			bool dominated = false;

			auto i = labels[v].begin();
			while(i != labels[v].end()) {

				Label * target_label = *i;

				if(comp_leq(target_label->point, new_cost)) {
					dominated = true;
					break;
				}

				if(comp_leq(new_cost, target_label->point)) {
					if(labels_in_queue.count(target_label) > 0) {
						target_label->mark_dominated = true;
						++i;
					} else {
						i = labels[v].erase(i);
						delete target_label;
					}
				} else {
					nondominatedLabels.emplace_back(target_label);
					++i;
				}
			}

			if(dominated) {
				delete new_cost;
				continue;
			}

			Label * new_label = new Label(new_cost, v, &*label);
			labels[v] = nondominatedLabels;
			labels[v].push_back(new_label);

			if(v != instance().target())
				lex_min_label.push(new_label);
		}

	}

	list<const Point *> target_labels;
	for(auto &label : labels[instance().target()])
		target_labels.push_back(new Point(*label->point));

	reset_solutions();
	add_solutions(target_labels.begin(), target_labels.end());

//	for(Label *label : labels[instance().target()]) {
//		const Label *current_label = label;
//		cout << "(";
//		while(current_label->n != instance().source()) {
//			cout << current_label->n << ", ";
//			current_label = current_label->pred;
//		}
//		cout << current_label->n << ")" << endl;
//	}

	node n;
	forall_nodes(n, instance().graph()) {
		for(auto &label : labels[n])
			delete label;
	}
}

}
