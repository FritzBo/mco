#pragma once
/*
 * est_dual_benson_scalarizer.h
 *
 *  Created on: 27.09.2013
 *      Author: fritz
 */

#ifndef EST_DUAL_BENSON_SCALARIZER_H_
#define EST_DUAL_BENSON_SCALARIZER_H_

#include <vector>
#include <list>

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>
#include <mco/basic/abstract_solver.h>
#include <mco/basic/abstract_graph_instance.h>
#include <mco/basic/weight_function_adaptors.h>
#include <mco/est/basic/kruskal_st_solver.h>
#include <mco/generic/benson_dual/dual_benson_scalarizer.h>
#include <mco/generic/benson_dual/ove_fp_v2.h>
#include <mco/generic/benson_dual/upper_image_container.h>

namespace mco {

class LexKruskalAdaptor {

public:

    LexKruskalAdaptor(AbstractGraphInstance& instance)
    :   instance_(instance) { }

    inline double operator()(const Point &weights, Point &value);

private:

    AbstractGraphInstance& instance_;

};

template<class OnlineVertexEnumerator = GraphlessOVE>
class ESTDualBensonScalarizer :

    public AbstractSolver<std::list<ogdf::edge>>,
    public AbstractUpperImageContainer<std::list<Point*>::const_iterator>

{

public:
	ESTDualBensonScalarizer(double epsilon = 1E-8)
    :   epsilon_(epsilon),
        cycles_(0) { }

	virtual ~ESTDualBensonScalarizer() {}

	void Solve(AbstractGraphInstance &graph) {

        for(auto point : inequalities_) {
            delete point;
        }

        for(auto point : extreme_points_) {
            delete point;
        }
        inequalities_.clear();
        extreme_points_.clear();

//        Point min_costs(numeric_limits<double>::infinity(),
//                        graph.dimension());
//
//        for(auto e : graph.graph().edges) {
//            for(unsigned i = 0; i < graph.dimension(); ++i) {
//                min_costs[i] = min(min_costs[i],
//                                   graph.weights()[e]->operator[](i));
//
//                if(min_costs[i] < 1) {
//                    min_costs[i] = 1;
//                }
//            }
//        }
//
//        ogdf::EdgeArray<Point*> normalized_edge_weights(graph.graph());
//
//        for(auto e : graph.graph().edges) {
//            auto p = new Point(graph.dimension());
//            for(unsigned i = 0; i < graph.dimension(); ++i) {
//                p->operator[](i) = graph.weights()[e]->operator[](i) / min_costs[i];
//            }
//            normalized_edge_weights[e] = p;
//        }

//        AbstractGraphInstance normalized_instance(graph.graph(),
//                                                  normalized_edge_weights,
//                                                  graph.dimension());

        std::list<Point *> frontier;

        instance_ = &graph;

        benson_scalarizer_ = new DualBensonScalarizer<OnlineVertexEnumerator>(LexKruskalAdaptor(graph),
                                    graph.dimension(),
                                    epsilon_);

//        benson_scalarizer_ = new DualBensonScalarizer<OnlineVertexEnumerator>(LexKruskalAdaptor(normalized_instance),
//                                                                              graph.dimension(),
//                                                                              epsilon_);

		benson_scalarizer_->Calculate_solutions(frontier);

        std::list<std::pair<std::list<ogdf::edge>, Point>> solutions;

        for(auto point : frontier) {
            solutions.push_back(make_pair(std::list<ogdf::edge>(), *point));
            extreme_points_.push_back(new Point(*point));
        }

        for(auto inequality : benson_scalarizer_->facets_list()) {
            inequalities_.push_back(new Point(inequality));
        }

        add_solutions(solutions.begin(), solutions.end());

//        for(auto e : graph.graph().edges) {
//            delete normalized_edge_weights[e];
//        }

	}

	double solver_time() {
		return cycles_ / (double) CLOCKS_PER_SEC;
	}

	double vertex_enumeration_time() {
		return benson_scalarizer_->vertex_enumeration_time();
	}

	int number_vertices() {
		return benson_scalarizer_->number_vertices();
	}

	int number_facets() {
		return benson_scalarizer_->number_facets();
	}

    virtual std::list<Point*>::const_iterator cbeginExtremePoints() {
        return extreme_points_.cbegin();
    }

    virtual std::list<Point*>::const_iterator cendExtremePoints() {
        return extreme_points_.cend();
    }

    virtual std::list<Point*>::const_iterator cbeginInequalities() {
        return inequalities_.cbegin();
    }

    virtual std::list<Point*>::const_iterator cendInequalities() {
        return inequalities_.cend();
    }

private:
	double epsilon_;
	DualBensonScalarizer<OnlineVertexEnumerator> *benson_scalarizer_;

	ogdf::NodeArray<ogdf::node> parents_;

	double kruskal_solver(std::vector<ogdf::edge> &sorted_edges, ogdf::EdgeArray<double> costs, Point &value);
	void make_set(ogdf::node);
	ogdf::node find_set(ogdf::node);
	void set_union(ogdf::node, ogdf::node);

	clock_t cycles_;

    AbstractGraphInstance* instance_;

    std::list<Point*> extreme_points_;
    std::list<Point*> inequalities_;

};

double LexKruskalAdaptor::
operator()(const Point &weights, Point &value) {
    LexKruskalSTSolver solver;

    solver.Solve(instance_.graph(),
                 LexWeightFunctionAdaptor(instance_.graph(),
                                          instance_.weights(),
                                          weights),
                 instance_.dimension() + 1);

    for(unsigned i = 0; i < instance_.dimension(); ++i) {
        value[i] = solver.min_cost()[i + 1];
    }

    return weights * value;
}

} /* namespace mco */
#endif /* EST_DUAL_BENSON_SCALARIZER_H_ */
