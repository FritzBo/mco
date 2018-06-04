/*
 * ep_solver_tsaggouris_approx.cpp
 *
 *  Created on: 26.03.2013
 *      Author: fritz
 */

#include <mco/ep/tsaggouris/ep_solver_tsaggouris_approx.h>

#include <vector>
#include <list>
#include <cmath>
#include <stdexcept>
#include <functional>

using std::invalid_argument;
using std::vector;
using std::list;
using std::log;
using std::floor;
using std::pair;
using std::max;
using std::min;
using std::numeric_limits;
using std::function;

#include <ogdf/basic/Graph.h>

using ogdf::node;
using ogdf::edge;
using ogdf::NodeArray;
using ogdf::EdgeArray;
using ogdf::Graph;

#include <mco/ep/basic/ep_instance.h>
#include <mco/ep/martins/label.h>
#include <mco/basic/forward_star.h>

namespace mco
{

unsigned EpSolverTsaggourisApprox::position(Instance &instance,
                                            const Point &epsilon,
                                            const Point &point) const
{

    unsigned int pos = 0;

    for(unsigned int i = 0; i < instance.dimension - 1; ++i)
    {
        double dpos = log(point[i] / c_min_[i]) / log(epsilon[i]);
        pos += bases_[i] * static_cast<unsigned>(floor(dpos));
    }
    return pos;
}

void EpSolverTsaggourisApprox::ExtendAndMerge(Instance &instance,
                                              vector<const Label *> &old_Py_n,
                                              const edge e,
                                              const node n,
                                              vector<const Label *> &new_Py_n,
                                              const Point &epsilon) const
{

    auto &weights = *instance.costs;
    auto dimension = instance.dimension;

    for(auto label : old_Py_n)
    {
        if(label == nullptr)
        {
            continue;
        }

        auto new_label = new const Label(label->point + *weights[e], n, label);

//        bool too_large = false;
//        for(unsigned i = 0; i < dimension - 1; ++i)
//        {
//            if(new_label->point[i] > max_path_cost_[i])
//            {
//                too_large = true;
//            }
//        }
//        if(too_large)
//        {
//            std::cout << new_label->point << std::endl;
//            const Label* curr = new_label;
//            unsigned length = 0;
//            while(curr->n->index() != 1 && length < 80)
//            {
//                std::cout << curr->n << ", ";
//                curr = curr->pred;
//                ++length;
//            }
//            std::cout << std::endl;
//            std::cout << "length: " << length << std::endl;
//            delete new_label;
//            continue;
//        }

        auto old_label = new_Py_n[position(instance, epsilon, new_label->point)];

        if(old_label == nullptr ||
           new_label->point[dimension - 1] < old_label->point[dimension - 1])
        {
            delete old_label;
            new_Py_n[position(instance, epsilon, new_label->point)] = new_label;
        }
        else
        {
            delete new_label;
        }
    }
}

void EpSolverTsaggourisApprox::Solve(Instance &instance,
                                     Point &epsilon)
{
    const unsigned int dimension = instance.dimension;
    const ForwardStar &graph = *instance.graph;
    const FSEdgeArray<Point *> &weights = *instance.costs;

    c_max_.resize(dimension - 1);
    c_min_.resize(dimension - 1);
    max_path_cost_.resize(dimension - 1);

    for(unsigned int i = 0; i < dimension - 1; ++i)
    {
        c_max_[i] = 0;
        c_min_[i] = numeric_limits<double>::infinity();
    }

    for(auto e : graph.edges)
    {
        for(unsigned int i = 0; i < dimension - 1; ++i)
        {
            c_max_[i] = std::max(c_max_[i], (*weights[e])[i]);
            c_min_[i] = std::min(c_min_[i], (*weights[e])[i]);

            if(c_min_[i] == 0)
            {
                throw invalid_argument("Edge cost 0 is not allowed.");
            }
        }
    }

    for(unsigned i = 0; i < dimension - 1; ++i)
    {
        max_path_cost_[i] = static_cast<unsigned>(c_max_[i] * (graph.numberOfNodes() - 1));
    }

//    std::cout << "c_max_: ";
//    for(auto x : c_max_)
//    {
//        std::cout << x << " ";
//    }
//    std::cout << std::endl;
//
//    std::cout << "c_min: ";
//    for(auto x : c_min_)
//    {
//        std::cout << x << " ";
//    }
//    std::cout << std::endl;

    double a;
    unsigned size = 1;
    bases_.resize(dimension - 1);

    for(unsigned i = 0; i < dimension - 1; ++i)
    {
        a = log(graph.numberOfNodes() * c_max_[i]) / log(epsilon[i]) + 1;
        bases_[i] = size;
        size *= static_cast<unsigned int>(ceil(a));
    }

//    std::cout << "bases: ";
//    for(auto x : bases_)
//    {
//        std::cout << x << " ";
//    }
//    std::cout << std::endl;

    auto *old_Py = new FSNodeArray<vector<const Label *>>(graph);

    for(auto n : graph.nodes)
    {
        (*old_Py)[n].resize(size, nullptr);
    }

    (*old_Py)[instance.source][0] = new Label(Point(dimension), -1, nullptr);

    for(unsigned i = 0; i < graph.numberOfNodes(); ++i)
    {
//        std::cout << i << std::endl;

        auto new_Py = new FSNodeArray<vector<const Label *>>(graph);
        for(auto n : graph.nodes)
        {
            (*new_Py)[n].resize(size, nullptr);
        }


        for(auto e : graph.edges)
        {
            ExtendAndMerge(instance,
                           (*old_Py)[graph.tail(e)],
                           e,
                           graph.tail(e),
                           (*new_Py)[graph.head(e)],
                           epsilon);
        }

        for(node n : graph.nodes)
        {
            for(auto l : old_Py->operator[](n))
            {
                delete l;
            }
        }
        delete old_Py;
        old_Py = new_Py;
    }

    vector<pair<csolution_type, const Point>> target_points;
    for(auto label : (*old_Py)[instance.target])
    {
        if(label != nullptr)
        {
            target_points.emplace_back(make_pair(list<edge>(), label->point));
        }
    }

    add_solutions(target_points.begin(), target_points.end());

    for(node n : graph.nodes)
    {
        for(auto l : old_Py->operator[](n))
        {
            delete l;
        }
    }
    delete old_Py;
}

} /* namespace mco */
