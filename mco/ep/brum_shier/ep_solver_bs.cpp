/*
 * bsssa.cpp
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#include <mco/ep/brum_shier/ep_solver_bs.h>

#include <queue>
#include <list>
#include <unordered_set>
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include <functional>

using std::queue;
using std::list;
using std::vector;
using std::unordered_set;
using std::cout;
using std::endl;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::node;
using ogdf::Graph;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using ogdf::AdjElement;

#include <mco/core/point.h>
#include <mco/ep/basic/ep_instance.h>

namespace mco {

/**
 * Returns the nondominated subset of the union of two sets of vectors of double arrays.
 * TODO: Using algorithm of Kung et al.
 */
bool DominationPartition(const list<const Point *> &source1,
                         const list<const Point *> &source2,
                         list<const Point *> &nondominated_subset,
                         list<const Point *> &dominated_subset,
                         double epsilon) {

//	unsigned int dim = source1.front()->dimension();

	bool new_labels = false;


//	if(dim == 2) {
//
//		auto it_s1 = source1.cbegin();
//		auto it_s2 = source2.cbegin();
//
//		if((**it_s1)[0] <= (**it_s2)[0] && (**it_s1)[1] > (**it_s2)[1]) {
//			nondominated_subset.push_back(*it_s1);
//			new_labels = true;
//			it_s1++;
//		} else {
//			nondominated_subset.push_back(*it_s2);
//			it_s2++;
//		}
//
//		while(!(it_s1 == source1.end() && it_s2 == source2.end())) {
//			if(it_s2 != source2.end() && (it_s1 == source1.end() || ((**it_s1)[0] >= (**it_s2)[0] && ((**it_s1)[1] < (**it_s2)[1])))) {
//				if((**it_s2)[1] <= (*nondominated_subset.back())[1] && (**it_s2) != (*nondominated_subset.back()))
//					nondominated_subset.push_back(*it_s2);
//				else
//					dominated_subset.push_back(*it_s2);
//				it_s2++;
//			} else {
//				if((**it_s1)[1] <= (*nondominated_subset.back())[1] && (**it_s1) != (*nondominated_subset.back())) {
//					nondominated_subset.push_back(*it_s1);
//					new_labels = true;
//				} else
//					dominated_subset.push_back(*it_s1);
//				it_s1++;
//			}
//		}
//
//	} else {

		vector<bool> marker(source2.size());
    
        EqualityPointComparator eq_comp(epsilon);
        ComponentwisePointComparator comp_leq(epsilon, false);
        ParetoDominationPointComparator pareto_dominates(epsilon);

		for(unsigned int i = 0; i < source2.size(); ++i)
			marker[i] = false;

		bool dominated;

		for(auto label_source1 : source1) {
			dominated = false;

			unsigned int i = 0;
			for(auto label_source2: source2) {

				if(eq_comp(label_source1, label_source2)) {
					dominated = true;
					break;
				}

				if(comp_leq(label_source1, label_source2)) {
					marker[i] = true;
                }

				if(comp_leq(label_source2, label_source2)) {
					dominated = true;
					break;
				}

				i++;
			}

			if(!dominated) {
				nondominated_subset.push_back(label_source1);
				new_labels = true;
			} else
				dominated_subset.push_back(label_source1);
		}

		unsigned int i = 0;
		for(auto label: source2) {
			if(!marker[i])
				nondominated_subset.push_back(label);
			else
				dominated_subset.push_back(label);
			i++;
		}
//	}

	for(auto dominated_label: dominated_subset) {
        
		dominated = false;
        
		for(auto nondominated_label: nondominated_subset)
            
			if(pareto_dominates(nondominated_label, dominated_label) ||
               eq_comp(nondominated_label, dominated_label)) {
                
				dominated = true;
				break;
			}
        
		if(!dominated) {
			std::cout << "Source1:" << std::endl;
			for(auto label: source1)
				std::cout << "(" << (*label)[0] << ", " << (*label)[1] << ")" << std::endl;

			std::cout << "Source2:" << std::endl;
			for(auto label: source2)
				std::cout << "(" << (*label)[0] << ", " << (*label)[1] << ")" << std::endl;

			std::cout << "Dominated:" << std::endl;
			for(auto label: dominated_subset)
				std::cout << "(" << (*label)[0] << ", " << (*label)[1] << ")" << std::endl;

			std::cout << "Nondominated:" << std::endl;
			for(auto label: nondominated_subset)
				std::cout << "(" << (*label)[0] << ", " << (*label)[1] << ")" << std::endl;

			std::cout << "New labels: " << new_labels << std::endl;

			assert(false);
		}
	}


	return new_labels;
}

void EpSolverBS::Solve() {
	const Graph &graph = instance().graph();
	const function<Point *(edge)> & weights(instance().weights());
	unsigned int dim = instance().dimension();
	const node source = instance().source();
	const node target = instance().target();

	queue<node> queue;
	unordered_set<node> nodes_in_queue;
	NodeArray<list<const Point *>> labels(graph);

	queue.push(source);
	nodes_in_queue.insert(source);

	labels[source].push_back(Point::Null(dim));

	while(!queue.empty()) {
		node n = queue.front();

//		cout << endl << n << ": ";

		list<const Point *> &currentNodeLabels = labels[n];

		AdjElement *adj;
		forall_adj(adj, n) {
			edge e = adj->theEdge();

			if(e->isSelfLoop())
				continue;

			node v = e->target();

			if(v == n)
				continue;

			if(v == source)
				continue;

//			cout << v << ", ";

			list<const Point *> new_labels;

			for(auto &label : currentNodeLabels) {
				Point * new_label = new Point(*label + *weights(e));
				new_labels.push_back(new_label);
			}

			if(labels[v].empty()) {

				labels[v].insert(labels[v].begin(), new_labels.begin(), new_labels.end());

				if(nodes_in_queue.count(v) == 0 && v != target) {
					queue.push(v);
					nodes_in_queue.insert(v);
				}

			} else {

				list<const Point *> dominated_subset;
				list<const Point *> nondominated_subset;

				bool changed = DominationPartition(new_labels, labels[v], nondominated_subset, dominated_subset, epsilon_);

				labels[v] = nondominated_subset;

				for(auto label : dominated_subset)
					delete label;

				if(changed && nodes_in_queue.count(v) == 0 && v != target) {
					queue.push(v);
					nodes_in_queue.insert(v);
				}
			}

		}

		queue.pop();
		nodes_in_queue.erase(n);

		assert(queue.size() <= static_cast<unsigned>(graph.numberOfNodes()));
	}

	node n;
	forall_nodes(n, graph)
		if(n != target)
			for(auto label : labels[n])
				delete label;

	add_solutions(labels[target].begin(), labels[target].end());
}

}
