//
//  label_selection.hpp
//  mco
//
//  Created by Fritz Bökler on 16.09.16.
//
//

#ifndef label_selection_h
#define label_selection_h

#include <list>
#include <functional>
#include <algorithm>

#include <mco/basic/abstract_solver.h>
#include <mco/basic/point.h>
#include <mco/basic/forward_star.h>

namespace mco {

class MOSPLabelSetting : public AbstractSolver<std::list<mco::node>>
{

public:
    struct InputDefinition {
        ForwardStar* graph;
        std::function<const Point*(const edge)> costs;
        unsigned dimension;
        node source;
        node target;
    };

    void Solve(InputDefinition& def);

private:

    struct NodeEntry;

    struct Label
    {
        unsigned label_id;
        Point cost;
        double sum;
        node n;
        edge pred_edge;
        unsigned pred_label_id;
        bool deleted;

        Label(const Label& label) = delete;
        void operator=(const Label& label) = delete;

        Label(Point&& cost,
              node n,
              edge p_edge,
              unsigned p_label_id,
              NodeEntry& node_entry)
        :   cost(std::move(cost)),
            n(n),
            pred_edge(p_edge),
            pred_label_id(p_label_id),
            deleted(false)
        {
            sum = 0;
            for(auto d : this->cost)
            {
                sum += d;
            }

            label_id = node_entry.get_new_id();
        }

        Label(Label&& other) noexcept
        :   label_id(other.label_id),
            cost(std::move(other.cost)),
            sum(other.sum),
            n(other.n),
            pred_edge(other.pred_edge),
            pred_label_id(other.pred_label_id),
            deleted(other.deleted)
        {
        }

        Label& operator=(Label&& other) noexcept
        {
            label_id = other.label_id;
            cost = std::move(other.cost);
            sum = other.sum;
            n = other.n;
            pred_edge = other.pred_edge;
            pred_label_id = other.pred_label_id;
            deleted = other.deleted;
            return *this;
        }

    };

    struct LabelComparator
    {
        bool operator()(const Label& l1, const Label& l2)
        {
            return l1.sum > l2.sum;
        }
    };

    struct NodeEntry
    {
        bool in_queue;

        inline NodeEntry()
        :   in_queue(false),
            closed_labels_(),
            open_labels_()
        {
        }

        NodeEntry(const NodeEntry&) = delete;
        NodeEntry(NodeEntry&&) = delete;

        NodeEntry& operator=(const NodeEntry&) = delete;
        NodeEntry& operator=(NodeEntry&&) = delete;

        inline void erase_closed(unsigned index)
        {
            closed_labels_[index] = std::move(closed_labels_[closed_labels_.size() - 1]);
            closed_labels_.pop_back();
        }

        inline void close_label(Label&& label)
        {
            closed_labels_.push_back(std::move(label));
        }

        inline void erase_open(unsigned index)
        {
            open_labels_[index].deleted = true;
        }

        inline void push_open(Label&& label)
        {
            open_labels_.push_back(std::move(label));
            std::push_heap(open_labels_.begin(),
                           open_labels_.end(),
                           LabelComparator());
        }

        inline Label& top_open()
        {
            return open_labels_.front();
        }

        inline void pop_open()
        {
            std::pop_heap(open_labels_.begin(),
                          open_labels_.end(),
                          LabelComparator());
            open_labels_.pop_back();
        }

        inline bool has_open()
        {
            return !open_labels_.empty();
        }

        inline std::vector<Label>& closed_labels()
        {
            return closed_labels_;
        }

        inline std::vector<Label>& open_labels()
        {
            return open_labels_;
        }

        unsigned get_new_id()
        {
            return id_counter++;
        }

    private:
        std::vector<Label> closed_labels_;
        std::vector<Label> open_labels_;

        unsigned id_counter = 0;
    };

    struct NodeComparator
    {
        NodeComparator(FSNodeArray<NodeEntry>& node_entries)
        :   node_entries_(node_entries) {}

        bool operator()(node n1, node n2)
        {
            return node_entries_[n1].top_open().sum >
                    node_entries_[n2].top_open().sum;
        }
    private:
        FSNodeArray<NodeEntry>& node_entries_;
    };

    bool non_dominated(Point new_cost,
                       NodeEntry& target_node_entry);
};

}

#endif /* label_selection_hpp */
