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
//#define STATS

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
#if defined USE_TREE_DELETION && defined STATS && !defined NDEBUG
        std::cout << "Tree-deleted labels: " << deleted_tree_labels_ << std::endl;
        std::cout << "Recursive deletions: " << recursive_deletions_ << std::endl;
#endif
#if defined STATS && !defined NDEBUG
        std::cout << "Processed recursively deleted labels: " << touched_recursively_deleted_label_ << std::endl;
        std::cout << "Number of compared labels: "  << label_compares_ / 1000 << "k" << std::endl;
        std::cout << "Numer of arc pushes: " << arc_pushes_ << std::endl;
#endif
    }

    unsigned deleted_tree_labels()
    {
        return deleted_tree_labels_;
    }

    unsigned recursive_deletions()
    {
        return recursive_deletions_;
    }

    unsigned touched_recursively_deleted_label()
    {
        return touched_recursively_deleted_label_;
    }

    unsigned long label_compares()
    {
        return label_compares_;
    }

    unsigned arc_pushes()
    {
        return arc_pushes_;
    }

    unsigned long deleted_labels()
    {
        return deleted_labels_;
    }

private:
    const double epsilon_;
    unsigned dimension_;

    bool            use_heuristic_  = false;
    heuristic_type  heuristic_;
    bool            use_bounds_     = false;
    Point           bounds_;

    unsigned long   label_compares_                     = 0;
    unsigned        deleted_tree_labels_                = 0;
    unsigned        recursive_deletions_                = 0;
    unsigned        arc_pushes_                         = 0;
    unsigned        touched_recursively_deleted_label_  = 0;
    unsigned long   deleted_labels_                     = 0;

    struct Label
    {
        Point cost;
        const ogdf::node n;
        const ogdf::edge pred_edge;
        Label* pred_label;
        bool deleted = false;
        bool in_queue = false;
        bool in_node_list = false;

#ifdef STATS
        bool mark_recursive_deleted = false;
#endif

#if defined USE_TREE_DELETION || defined STATS
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

#if defined USE_TREE_DELETION || defined STATS
    void recursive_delete(Label& label);
#endif

    bool check_heuristic_prunable(const Label& label);


};

}

#endif /* __mco__ep_lc_ls_ */
