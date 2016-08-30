/*
 * bsssa.cpp
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#include <mco/ep/brum_shier/ep_solver_bs_bi.h>

#include <list>
#include <deque>
#include <iostream>
#include <memory>
#include <vector>
#include <cassert>
#include <functional>

using std::list;
using std::vector;
using std::cout;
using std::endl;
using std::function;
using std::pair;

#include <mco/basic/point.h>
#include <mco/ep/basic/ep_instance.h>
#include <mco/basic/utility.h>
#include <mco/basic/ringbuffer.h>

namespace mco {

bool EpSolverBSBi::
check_domination(vector<Label>& new_labels,
                 NodeEntry& neighbor_entry) {

    EqualityPointComparator eq;
    ComponentwisePointComparator leq(0, false);

    auto& neighbor_labels = neighbor_entry.labels();
    auto& neighbor_new_labels = neighbor_entry.new_labels();

//    bool changed = false;

//    merge(new_labels, neighbor_labels);
//    changed = merge(new_labels, neighbor_new_labels, true);

    merge(new_labels, neighbor_labels, neighbor_new_labels);

    assert(is_sorted(new_labels.cbegin(),
                     new_labels.cend()));
    assert(is_sorted(neighbor_labels.cbegin(),
                     neighbor_labels.cend()));
    assert(is_sorted(neighbor_new_labels.cbegin(),
                     neighbor_new_labels.cend()));

    return !neighbor_new_labels.empty();
}

#if defined USE_TREE_DELETION || defined STATS
void EpSolverBSBi::
recursive_delete(Label& label)
{
#if defined STATS && defined USE_TREE_DELETION
    recursive_deletions_ += 1;
#endif

    std::deque<Label*> queue;
    queue.push_back(&label);

    while(!queue.empty()) {
        Label* curr = queue.front();
        queue.pop_front();

#ifdef USE_TREE_DELETION
        curr->deleted = true;
#else
        curr->mark_recursive_deleted = true;
#endif

#if defined USE_TREE_DELETION && defined STATS
        ++deleted_tree_labels_;
#endif

        for(auto succ : curr->successors) {
            queue.push_back(succ);
        }
    }
}
#endif

bool EpSolverBSBi::check_heuristic_prunable(const Label& label) {

    if(label.cost.first + heuristic_(label.n, 0) > bounds_[0])
    {
        return true;
    }
    else if(label.cost.second + heuristic_(label.n, 1) > bounds_[1])
    {
        return true;
    }

    return false;
}


void EpSolverBSBi::Solve(const ForwardStar& graph,
                       std::function<const BiPoint*(const edge)> weights,
                       unsigned dimension,
                       const node source,
                       const node target) {

    assert(dimension == 2);

    dimension_ = dimension;

    ring_buffer<node> queue(graph.numberOfNodes() + 1);
	FSNodeArray<NodeEntry> node_entries(graph);

    node_entries[source].new_labels().emplace_back(BiPoint(0, 0),
                                                   source,
                                                   nulledge,
                                                   0,
                                                   node_entries[source]);

	queue.push_back(source);
    node_entries[source].in_queue = true;

	while(!queue.empty()) {
		node current_node = queue.front();
        queue.pop_front();

//		cout << n << ": ";

		NodeEntry& current_node_entry = node_entries[current_node];
        current_node_entry.in_queue = false;

        if(current_node_entry.has_new_labels()) {

            auto& current_new_labels = current_node_entry.new_labels();

            vector<Label> new_labels;

            for(auto current_edge : graph.adj_edges(current_node))
            {
                node neighbor = graph.head(current_edge);

                if(neighbor == source)
                {
                    continue;
                }

#ifdef STATS
                arc_pushes_ += 1;
#endif

                auto& neighbor_entry = node_entries[neighbor];

    //			cout << neighbor << ", ";

                new_labels.resize(0);

                size_t size = 0;

                for(auto& label : current_new_labels)
                {
#ifdef STATS
                    if(label.mark_recursive_deleted)
                    {
                        ++touched_recursively_deleted_label_;
                    }
#endif

                    assert(label.n == current_node);

#if defined USE_TREE_DELETION || defined STATS
                    if(!label.deleted)
                    {
#endif

                        Label new_label(label.cost + *weights(current_edge),
                                        neighbor,
                                        current_edge,
                                        label.label_id,
                                        neighbor_entry);

                        if(!use_bounds_ || !use_heuristic_ ||
                           !check_heuristic_prunable(new_label))
                        {

                            new_labels.push_back(std::move(new_label));
                            ++size;
                        }
#ifdef STATS
                        else
                        {
                            ++deleted_labels_;
                        }
#endif
#if defined USE_TREE_DELETION || defined STATS
                    }
#endif
                }

                if(!new_labels.empty())
                {
                    bool changed = check_domination(new_labels,
                                                    neighbor_entry);

                    if(!node_entries[neighbor].in_queue &&
                       changed &&
                       neighbor != target &&
                       neighbor != source) {

                        queue.push_back(neighbor);
                        node_entries[neighbor].in_queue = true;
                    }

                }
            }

            current_node_entry.proceed_labels_it();

        }

    }

    auto& target_entry = node_entries[target];
    for(size_t i = 0; i < target_entry.new_labels().size(); ++i)
    {
        auto& label = target_entry.new_labels()[i];
        list<node> path;
        auto curr = &label;

#if defined USE_TREE_DELETION || defined STATS
        if(!curr->deleted)
        {
#endif
            while(curr->n != source)
            {
                path.push_back(curr->n);

#ifndef NDEBUG
                Label* test = curr;
#endif
                node pred_node = graph.tail(curr->pred_edge);

                assert(pred_node != curr->n);

                for(auto& candidate_label : node_entries[pred_node].labels())
                {
                    if(candidate_label.label_id == curr->pred_label_id)
                    {
                        curr = &candidate_label;
                        break;
                    }
                }

                assert(test != curr);
            }

            assert(path.size() > 0);

            path.push_back(source);

            path.reverse();

//            add_solution(path, label.cost);
            add_solution(path, Point({static_cast<double>(label.cost.first),
                                        static_cast<double>(label.cost.second)}));
#if defined USE_TREE_DELETION || defined STATS
        }
#endif
    }
}

int EpSolverBSBi::compare(Label const &a, Label const &b) {
    if (a.cost.first <= b.cost.first && a.cost.second < b.cost.second) return -1;
    if (a.cost.first < b.cost.first && a.cost.second <= b.cost.second) return -1;
    if (a.cost.first >= b.cost.first && a.cost.second >= b.cost.second) return 1;
//    if (a.cost.first > b.cost.first && a.cost.second >= b.cost.second) return 1;
    return 0;
}

void EpSolverBSBi::merge(vector<Label> &N,
                         vector<Label> &pushed,
                         vector<Label> &pending) {

    bool print_pushed = false;

    vector<bool> nDeleted(N.size(), false);
    int nIndex = 0;
    auto nIt = N.begin();
    auto pushedIt = pushed.begin();
    int deletedCount = 0;
    while (nIt != N.end() && pushedIt != pushed.end()) {
        switch(compare(*nIt, *pushedIt)) {
            case -1:
                ++pushedIt;
                ++deletedCount;
                break;
            case 1:
                nDeleted[nIndex] = true;
                ++nIt;
                ++nIndex;
                break;
            case 0:
                if (nIt->cost.first > pushedIt->cost.first) {
                    using std::iter_swap;
                    iter_swap(pushedIt, pushedIt - deletedCount);
                    ++pushedIt;
                } else {
                    ++nIt;
                    ++nIndex;
                }
        }
    }

    while (pushedIt != pushed.end()) {
        using std::iter_swap;
        iter_swap(pushedIt, pushedIt - deletedCount);
        ++pushedIt;
    }
    pushed.resize(pushed.size() - deletedCount);

    if(print_pushed)
    {
        for(auto& l : pushed)
        {
            std::cout << l.cost.first << ", " << l.cost.second << std::endl;
        }
    }

    auto pendingIt = pending.begin();
    nIt = N.begin();
    nIndex = 0;
    while (nIt != N.end() && nDeleted[nIndex]) {
        ++nIt;
        ++nIndex;
    }
    deletedCount = 0;
    while (nIt != N.end() && pendingIt != pending.end()) {
        switch(compare(*nIt, *pendingIt)) {
            case -1:
                ++pendingIt;
                ++deletedCount;
                break;
            case 1:
                nDeleted[nIndex] = true;
                ++nIt;
                ++nIndex;
                while (nIt != N.end() && nDeleted[nIndex]) {
                    ++nIt;
                    ++nIndex;
                }
                break;
            case 0:
                if (nIt->cost.first > pendingIt->cost.first) {
                    using std::iter_swap;
                    iter_swap(pendingIt, pendingIt - deletedCount);
                    ++pendingIt;
                } else {
                    ++nIt;
                    ++nIndex;
                }
        }
    }

    while (pendingIt != pending.end()) {
        using std::iter_swap;
        iter_swap(pendingIt, pendingIt - deletedCount);
        ++pendingIt;
    }
    pending.resize(pending.size() - deletedCount);
    size_t oldPendingCount = pending.size();

    nIndex = 0;
    for (auto& a : N) {
        if (!nDeleted[nIndex]) {
            pending.push_back(std::move(a));
        }
        ++nIndex;
    }

    std::inplace_merge(pending.begin(), pending.begin() + oldPendingCount, pending.end(), [] (const Label& a, const Label& b) {return a.cost < b.cost;});
}


//bool EpSolverBSBi::merge(vector<Label>& new_set,
//                         vector<Label>& target_set,
//                         bool insert)
//{
//    if(target_set.empty())
//    {
//        if(!insert)
//        {
//            return false;
//        }
//        else
//        {
//            for(auto& l : new_set)
//            {
//                target_set.push_back(std::move(l));
//            }
//            return true;
//        }
//    }
//
//    bool changed = false;
//    unsigned new_idx     = 0;
//    unsigned target_idx  = 0;
//
//    vector<unsigned> new_set_deleted;
//    vector<unsigned> target_set_deleted;
//    vector<pair<unsigned, unsigned>> insertions;
//
//    while(new_idx < new_set.size() || target_idx < target_set.size())
//    {
//        if(new_idx == new_set.size())
//        {
//            break;
//        }
//        else if(target_idx == target_set.size())
//        {
//            if(insert)
//            {
//                auto& new_label = new_set[new_idx];
//                auto& last_target_label = target_set[target_set.size() - 1];
//
//                while(new_idx < new_set.size())
//                {
//                    if((new_label.cost[0] >= last_target_label.cost[0] &&
//                       new_label.cost[1] < last_target_label.cost[1]) ||
//                       (new_label.cost[0] > last_target_label.cost[0] &&
//                       new_label.cost[1] <= last_target_label.cost[1]))
//                    {
//                        target_set.push_back(std::move(new_set[new_idx]));
//                        changed = true;
//                    }
//                    ++new_idx;
//                }
//            }
//            break;
//        }
//
//        auto& new_label = new_set[new_idx];
//        auto& target_label = target_set[target_idx];
//
//        if(new_label.cost[0] <= target_label.cost[0] &&
//           new_label.cost[1] <= target_label.cost[1])
//        {
//            target_set_deleted.push_back(target_idx);
//            ++target_idx;
//        }
//        else if(new_label.cost[0] >= target_label.cost[0] &&
//                new_label.cost[1] >= target_label.cost[1])
//        {
//            new_set_deleted.push_back(new_idx);
//            ++new_idx;
//        }
//        else if(new_label.cost[0] < target_label.cost[0] &&
//                new_label.cost[1] >= target_label.cost[1])
//        {
//            if(insert)
//            {
//                insertions.push_back(std::make_pair(new_idx, target_idx));
//                changed = true;
//            }
//            ++new_idx;
//        }
//        else
//        {
//            ++target_idx;
//        }
//    }
//
//    if(!new_set_deleted.empty())
//    {
//        clean_set(new_set, new_set_deleted, &insertions);
//    }
//    
//    if(!target_set_deleted.empty())
//    {
//        clean_set(target_set, target_set_deleted, nullptr);
//    }
//
//    if(!insertions.empty())
//    {
//        unsigned shift = insertions.size();
//        target_set.resize(target_set.size() + insertions.size());
//
//        auto insertion_iter = insertions.rbegin();
//
//        for(unsigned i = target_set.size() - 1; i >= insertions[0].second; ++i)
//        {
//            if(insertion_iter->second == i)
//            {
//                target_set[i + shift] = std::move(new_set[insertion_iter->first]);
//                ++insertion_iter;
//                --shift;
//            }
//            else
//            {
//                target_set[i + shift] = std::move(target_set[i]);
//            }
//        }
//    }
//
//    return changed;
//}
//
//void EpSolverBSBi::clean_set(std::vector<Label>& set,
//                             std::vector<unsigned>& positions,
//                             std::vector<pair<unsigned, unsigned>>* insertions)
//{
//    unsigned shift = 0;
//    auto cur_position = positions.begin();
//
//    vector<pair<unsigned, unsigned>>::iterator cur_insertion;
//    if(insertions != nullptr)
//    {
//        cur_insertion = insertions->begin();
//    }
//
//    for(unsigned i = 0; i < set.size(); ++i)
//    {
//        if(insertions != nullptr && cur_insertion != insertions->end())
//        {
//            if(cur_insertion->first == i)
//            {
//                cur_insertion->first -= shift;
//                ++cur_insertion;
//            }
//        }
//        if(shift > 0)
//        {
//            set[i - shift] = std::move(set[i]);
//        }
//        if(cur_position != positions.end() && *cur_position == i)
//        {
//            ++shift;
//            ++cur_position;
//        }
//    }
//    for(unsigned i = 0; i < shift; ++i)
//    {
//        set.pop_back();
//    }
//}
}
