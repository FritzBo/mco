#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_BI_H_
#define BSSSA_BI_H_

//#define USE_TREE_DELETION
//#define STATS

#include <iostream>

#include <mco/basic/abstract_solver.h>

#include <mco/basic/forward_star.h>

namespace mco
{

using BiPoint = std::pair<int, int>;

class EpSolverBSBi : public AbstractSolver<std::list<mco::node>>
{

    using heuristic_type = std::function<double(node, unsigned)>;
    
public:
	EpSolverBSBi(double epsilon = 0)
    : epsilon_(epsilon)
    {
    }
    
	virtual void Solve(const ForwardStar& graph,
                       std::function<const BiPoint*(const edge)> costs,
                       unsigned dimension,
                       const node source,
                       const node target);

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

    ~EpSolverBSBi()
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

    bool use_heuristic_ = false;
    heuristic_type heuristic_;
    bool use_bounds_ = false;
    Point bounds_;

    unsigned deleted_tree_labels_ = 0;
    unsigned recursive_deletions_ = 0;

    unsigned touched_recursively_deleted_label_ = 0;
    unsigned long label_compares_ = 0;
    unsigned arc_pushes_ = 0;
    unsigned long deleted_labels_ = 0;

    struct NodeEntry;

    struct Label
    {
        unsigned label_id;
        BiPoint cost;
        node n;
        edge pred_edge;
        unsigned pred_label_id;

#ifdef STATS
        bool mark_recursive_deleted = false;
#endif

#if defined USE_TREE_DELETION || defined STATS
        /* Lists are more efficient here, because most of the time
         the collection of successors is actually empty or very small. */
        std::list<Label*> successors;
        bool deleted = false;
#endif

        Label()
        :   cost(BiPoint()),
            n(nullnode),
            pred_edge(nulledge),
            pred_label_id(0)
        {
            label_id = 0;
        }

        Label(BiPoint cost,
              node n,
              edge p_edge,
              unsigned p_label_id,
              NodeEntry& node_entry)
        :   cost(cost),
            n(n),
            pred_edge(p_edge),
            pred_label_id(p_label_id)

#if defined USE_TREE_DELETION || defined STATS
            , deleted(false)
#endif
        {
            label_id = node_entry.get_new_id();
        }

        Label(Label&& other)
        :   label_id(other.label_id),
            cost(std::move(other.cost)),
            n(other.n),
            pred_edge(other.pred_edge),
            pred_label_id(other.pred_label_id)
#if defined USE_TREE_DELETION || defined STATS
            , deleted(false)
            , successors(std::move(other.successors))
#endif
        {
        }

        Label& operator=(Label&& other)
        {
            label_id = other.label_id;
            cost = std::move(other.cost);
            n = other.n;
            pred_edge = other.pred_edge;
            pred_label_id = other.pred_label_id;
#if defined USE_TREE_DELETION || defined STATS
            deleted = other.deleted;
            successors = std::move(other.successors);
#endif
            return *this;
        }
    };

    struct NodeEntry
    {
        using value_type = Label;

        bool in_queue;

        inline NodeEntry()
        :   in_queue(false),
            labels_(),
            new_labels_()
        {
        }

        inline void erase(std::vector<Label>& arr,
                   unsigned index)
        {

            assert(&arr == &labels_ ||
                   &arr == &new_labels_);

            arr[index] = std::move(arr[arr.size() - 1]);
            arr.pop_back();
        }

        inline void proceed_labels_it()
        {
            for(auto& label : new_labels_)
            {
                labels_.push_back(std::move(label));
            }

            new_labels_.clear();
        }

        inline bool has_new_labels()
        {
            return !new_labels_.empty();
        }

        inline std::vector<Label>& labels()
        {
            return labels_;
        }

        inline std::vector<Label>& new_labels()
        {
            return new_labels_;
        }

        unsigned get_new_id()
        {
            return id_counter++;
        }

    private:
        std::vector<Label> labels_;
        std::vector<Label> new_labels_;

        unsigned id_counter = 1;
    };

    bool check_domination(std::vector<Label>& new_labels,
                          NodeEntry& neighbor_entry);

#if defined USE_TREE_DELETION || defined STATS
    void recursive_delete(Label& label);
#endif

    bool check_heuristic_prunable(const Label& label);

    inline void remove_label(Label& label);

    int compare(Label const &a, Label const &b);

    void merge(std::vector<Label> &N,
               std::vector<Label> &pushed,
               std::vector<Label> &pending);

//    bool merge(std::vector<Label>& input_set,
//               std::vector<Label>& target_set,
//               bool insert = false);
//
//    void clean_set(std::vector<Label>& set,
//                   std::vector<unsigned>& positions,
//                   std::vector<std::pair<unsigned, unsigned>>* insertions);

};

inline void EpSolverBSBi::remove_label(Label& label)
{
#if defined USE_TREE_DELETION || defined STATS
#if defined STATS && !defined USE_TREE_DELETION
    if(label.pred_label != nullptr)
    {
#endif

        auto it = find(label.pred_label->successors.begin(),
                       label.pred_label->successors.end(),
                       label);

        assert(it != label.pred_label->successors.end());

        label.pred_label->successors.erase(it);

#if defined STATS && !defined USE_TREE_DELETION
    }

    for(auto label : label->successors)
    {
        label->pred_label = nullptr;
    }
#endif


    if(label.successors.empty())
    {
        label.deleted = true;
    }
    else
    {
        recursive_delete(*label);
    }
#endif
}

inline BiPoint operator+(const BiPoint& p1, const BiPoint& p2)
{
    BiPoint p;
    p.first = p1.first + p2.first;
    p.second = p1.second + p2.second;
    return p;
}

}

#endif /* BSSSA_H_ */
