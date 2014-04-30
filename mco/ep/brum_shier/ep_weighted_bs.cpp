/*
 * bsssa.cpp
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#include <mco/ep/brum_shier/ep_weighted_bs.h>

#include <queue>
#include <list>
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include <functional>

using std::queue;
using std::list;
using std::vector;
using std::cout;
using std::endl;
using std::function;

#include <setoper.h>
#include <cdd.h>
#include <ogdf/basic/Graph.h>

using ogdf::edge;
using ogdf::node;
using ogdf::Graph;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using ogdf::AdjElement;

#include <mco/basic/point.h>
#include <mco/ep/basic/ep_instance.h>
#include <mco/basic/utility.h>

namespace mco {

bool EpWeightedBS::
ConvexHull(const list<const Point *> &source1,
           const list<const Point *> &source2,
           list<const Point *> &nondominated_subset,
           list<const Point *> &dominated_subset,
           double epsilon) {
    
    unsigned dim = source1.back()->dimension();
    unsigned no_points = source1.size() + source2.size();
    bool changed = false;
    
    dd_MatrixPtr points = dd_CreateMatrix(no_points + dim,
                                          dim + 1);
    
    points->representation = dd_Generator;
    
    int row_index = 0;
    
    // Insert extreme rays
    for(unsigned j = 0; j < dim; ++j) {
        
        for(unsigned k = 1; k < dim + 1; ++k) {
            dd_set_d(points->matrix[row_index][k], j + 1 == k ? 1 : 0);
        }
        
        dd_set_d(points->matrix[row_index][0], 0);
        
        ++row_index;
    }
    
    // Insert points from source1 (new points)
    for(auto pointPtr : source1) {
        Point point = *pointPtr;
        
        for(unsigned k = 1; k < dim + 1; ++k) {
            dd_set_d(points->matrix[row_index][k], point[k - 1]);
        }
        
        dd_set_d(points->matrix[row_index][0], 1);
        
        ++row_index;
                                    
    }
                                
    // Insert points from source2 (old points)
    for(auto pointPtr : source2) {
        Point point = *pointPtr;
        
        for(unsigned k = 1; k < dim + 1; ++k) {
            dd_set_d(points->matrix[row_index][k], point[k - 1]);
        }
        
        dd_set_d(points->matrix[row_index][0], 1);
        
        ++row_index;
        
    }
    
//    cout << no_points << endl;
    
//    dd_WriteMatrix(stdout, points);

    dd_ErrorType err;
    dd_rowset redundand_rows;
    
    redundand_rows = dd_RedundantRows(points, &err);
    
    auto iter = source1.begin();
    
    for(unsigned i = dim; i < row_index; ++i) {
        
        if(iter == source1.end()) {
            iter = source2.begin();
        }
        
        if(set_member(i + 1, redundand_rows)) {
            dominated_subset.push_back(*iter);
        } else {
            nondominated_subset.push_back(*iter);
        }
        
        ++iter;
    }

    dd_FreeMatrix(points);
    
    EqualityPointComparator eq;
    
    iter = nondominated_subset.begin();
    for(auto p : source2) {
        if(iter == nondominated_subset.end()) {
            break;
        }
        
        if(!eq(p, *iter)) {
            changed = true;
        }
        
        ++iter;
    }
    
    // FIXME
//    if(changed) {
//        cout << "old" << endl;
//        for(auto p: source2) {
//            cout << *p << endl;
//        }
//        cout << "new" << endl;
//        for(auto p: nondominated_subset) {
//            cout << *p << endl;
//        }
//        cout << endl;
//    }

    return changed;
}

void EpWeightedBS::Solve(const Graph& graph,
                       std::function<const Point*(const ogdf::edge)> weights,
                       unsigned dim,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed) {
    
	queue<node> queue;
	NodeArray<bool> nodes_in_queue(graph, false);
	NodeArray<list<const Point *>> labels(graph);
    
    dd_set_global_constants();

	queue.push(source);
	nodes_in_queue[source] = true;

	labels[source].push_back(Point::Null(dim));

	while(!queue.empty()) {
		node n = queue.front();

//		cout << n << ": ";

		list<const Point *> &currentNodeLabels = labels[n];

        for(auto adj : n->adjEdges) {
            
            edge e = adj->theEdge();
            
			if(e->isSelfLoop()) {
				continue;
            }
            
            node v;
            
            if(directed) {

                v = e->target();

                if(v == n)
                    continue;
                
            } else {
                
                v = e->target() == n ? e->source() : e->target();
            }
            
            if(v == source) {
                continue;
            }

//			cout << v << ", ";

			list<const Point *> new_labels;

			for(auto &label : currentNodeLabels) {
				Point * new_label = new Point(*label + *weights(e));
				new_labels.push_back(new_label);
			}

			if(labels[v].empty()) {

				labels[v].insert(labels[v].begin(), new_labels.begin(), new_labels.end());

				if(!nodes_in_queue[v] && v != target) {
					queue.push(v);
					nodes_in_queue[v] = true;
				}

			} else {

				list<const Point *> dominated_subset;
				list<const Point *> nondominated_subset;

				bool changed = ConvexHull(new_labels,
                                          labels[v],
                                          nondominated_subset,
                                          dominated_subset,
                                          epsilon_);

				labels[v] = nondominated_subset;

				for(auto label : dominated_subset) {
					delete label;
                }

				if(changed && !nodes_in_queue[v] && v != target) {
					queue.push(v);
					nodes_in_queue[v] = true;
				}
                
			}

		}

		queue.pop();
		nodes_in_queue[n] = false;
        
//        cout << endl;

		assert(queue.size() <= static_cast<unsigned>(graph.numberOfNodes()));
	}

    for(auto n : graph.nodes) {
		if(n != target) {
			for(auto label : labels[n]) {
				delete label;
            }
        }
    }

	add_solutions(labels[target].begin(), labels[target].end());
    
    dd_free_global_constants();
}

}
