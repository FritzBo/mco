#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_PC_H_
#define BSSSA_PC_H_

//#define USE_TREE_DELETION
//#define STATS

#include <iostream>

#include <mco/basic/abstract_solver.h>
#include <mco/basic/ringbuffer.h>
#include <mco/basic/forward_star.h>

namespace mco
{

class EpSolverBSPc : public AbstractSolver<std::list<mco::node>>
{

    using heuristic_type = std::function<double(node, unsigned)>;
    
public:
	EpSolverBSPc(double epsilon = 0)
    : epsilon_(epsilon)
    {
    }
    
	virtual void Solve(const ForwardStar& graph,
                       std::function<const Point*(const edge)> costs,
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

    ~EpSolverBSPc()
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
        Point cost;
        node n;
        edge pred_edge;
        Label* pred_label;

#ifdef STATS
        bool mark_recursive_deleted = false;
#endif

#if defined USE_TREE_DELETION || defined STATS
        /* Lists are more efficient here, because most of the time
         the collection of successors is actually empty or very small. */
        bool deleted = false;
#endif

        std::list<Label*> successors;

        Label(Point&& cost,
              node n,
              edge p_edge,
              Label* p_label)
        :   cost(std::move(cost)),
            n(n),
            pred_edge(p_edge),
            pred_label(p_label)
#if defined USE_TREE_DELETION || defined STATS
            , deleted(false)
#endif
        {
        }

        Label(Label&& other)
        :   cost(std::move(other.cost)),
            n(other.n),
            pred_edge(other.pred_edge),
            pred_label(other.pred_label)
#if defined USE_TREE_DELETION || defined STATS
            , deleted(false)
            , successors(std::move(other.successors))
#endif
        {
        }

        Label& operator=(Label&& other)
        {
            cost = std::move(other.cost);
            n = other.n;
            pred_edge = other.pred_edge;
            pred_label = other.pred_label;
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

        inline void erase(std::vector<Label*>& arr,
                   unsigned index)
        {

            assert(&arr == &labels_ ||
                   &arr == &new_labels_);

            auto erase_label = arr[index];

            for(auto succ_label : erase_label->successors)
            {
                succ_label->pred_label = nullptr;
            }

            if(erase_label->pred_label != nullptr)
            {
                auto it = erase_label->pred_label->successors.begin();
                while(it != erase_label->pred_label->successors.end())
                {
                    if(*it == erase_label)
                    {
                        erase_label->pred_label->successors.erase(it);
                        break;
                    }
                    ++it;
                }
            }

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

        inline std::vector<Label*>& labels()
        {
            return labels_;
        }

        inline std::vector<Label*>& new_labels()
        {
            return new_labels_;
        }

    private:
        std::vector<Label*> labels_;
        std::vector<Label*> new_labels_;
    };

    bool check_domination(std::vector<Label*>& new_labels,
                          NodeEntry& neighbor_entry);

#if defined USE_TREE_DELETION || defined STATS
    void recursive_delete(Label& label);
#endif

    bool check_heuristic_prunable(const Label& label);

    inline void remove_label(Label& label);

    bool parent_check(node current_node,
                      FSNodeArray<NodeEntry>& node_entries,
                      const ForwardStar& graph);
};

inline void EpSolverBSPc::remove_label(Label& label)
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
}

#endif /* BSSSA_PC_H_ */
