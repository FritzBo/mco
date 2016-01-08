#pragma once
/*
 * ep_lc_ls.h
 *
 *  Created on: 1.12.2015
 *      Author: fritz
 */

#ifndef __mco__ep_lc_ls_
#define __mco__ep_lc_ls_

//#define USE_TREE_DELETION
#define STATS

#include <mco/basic/abstract_solver.h>

namespace mco
{

class EpLcLs : public AbstractSolver<std::list<ogdf::edge>>
{

    using heuristic_type = std::function<double(ogdf::node, unsigned)>;
    
public:
	EpLcLs(double epsilon = 0)
    : epsilon_(epsilon) { }
    
	virtual void Solve(const ogdf::Graph& graph,
                       std::function<const Point*(const ogdf::edge)> costs,
                       unsigned dimension,
                       const ogdf::node source,
                       const ogdf::node target,
                       bool directed = true);

    void set_heuristic(heuristic_type heuristic)
    {
        heuristic_ = heuristic;
        use_heuristic_ = true;
    }

    void set_bounds(Point bounds)
    {
        bounds_ = std::move(bounds);
        use_bounds_ = true;
    }

    ~EpLcLs()
    {
#ifdef STATS
        std::cout << "Number of compared labels: " << label_compares_ / 1000 << "k" << std::endl;
#endif
    }


private:
    const double epsilon_;
    unsigned dimension_;

    bool use_heuristic_ = false;
    heuristic_type heuristic_;
    bool use_bounds_ = false;
    Point bounds_;

#ifdef STATS
    unsigned long label_compares_ = 0;
#endif

    struct Label
    {
        Point cost;
        const ogdf::node n;
        const ogdf::edge pred_edge;
        Label* pred_label;
        bool deleted = false;
        bool in_queue = false;
        bool in_node_list = false;

#ifdef USE_TREE_DELETION
        std::list<Label*> successors;
#endif

        Label(const Point cost,
              const ogdf::node n,
              const ogdf::edge p_edge,
              Label* p_label)
        :   cost(cost),
            n(n),
            pred_edge(p_edge),
            pred_label(p_label),
            deleted(false)
        {
        }

    };

    struct NodeEntry
    {

        NodeEntry()
        :   labels_(),
            labels_end_(0)
        {

        }

        inline void push_back(Label* label) {

            if(labels_.size() > labels_end_)
            {
                labels_[labels_end_] = label;
            }
            else
            {
                labels_.push_back(label);
            }

            ++labels_end_;

        }

        inline void erase(unsigned index)
        {

            if(labels_[index]->in_queue)
            {
                labels_[index]->deleted = true;
            }
            else
            {
                delete labels_[index];
            }

            labels_[index]->in_node_list = false;

            labels_[index] = labels_[labels_end_ - 1];

#ifndef NDEBUG
            labels_[labels_end_ - 1] = nullptr;
#endif
            --labels_end_;
        }

        inline std::vector<Label*>& labels()
        {
            return labels_;
        }

        inline unsigned labels_end()
        {
            return labels_end_;
        }

        ~NodeEntry()
        {
            for(unsigned i = 0; i < labels_end_; ++i)
            {
                delete labels_[i];
            }
        }

    private:
        std::vector<Label*> labels_;
        unsigned labels_end_;
    };

    bool check_domination(Label* new_label,
                          NodeEntry& neighbor_entry);

#ifdef USE_TREE_DELETION
    void recursive_delete(Label& label);
#endif

    bool check_heuristic_prunable(const Label& label);


};

}

#endif /* __mco__ep_lc_ls_ */
