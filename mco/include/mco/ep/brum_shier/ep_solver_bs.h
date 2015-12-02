#pragma once
/*
 * bsssa.h
 *
 *  Created on: 13.03.2013
 *      Author: fritz
 */

#ifndef BSSSA_H_
#define BSSSA_H_

#include <mco/basic/abstract_solver.h>

namespace mco
{

class EpSolverBS : public AbstractSolver<std::list<ogdf::edge>>
{

    using heuristic_type = std::function<double(ogdf::node, unsigned)>;
    
public:
	EpSolverBS(double epsilon = 0)
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


private:
    const double epsilon_;
    unsigned dimension_;

    bool use_heuristic_ = false;
    heuristic_type heuristic_;
    bool use_bounds_ = false;
    Point bounds_;

    struct Label
    {
        const Point cost;
        const ogdf::node n;
        const ogdf::edge pred_edge;
        Label* pred_label;
        bool deleted = false;

        std::list<Label*> successors;

        Label(const Point cost,
              const ogdf::node n,
              const ogdf::edge p_edge,
              Label* p_label)
        :   cost(cost),
            n(n),
            pred_edge(p_edge),
            pred_label(p_label),
            deleted(false) { }

//        void erase_successor(const Label* label)
//        {
//#ifndef NDEBUG
//            bool found = false;
//#endif
//            for(unsigned i = 0; i < successors_.size(); ++i)
//            {
//                if(successors_[i] == label)
//                {
//                    successors_[i] = successors_[successors_.size() - 1];
//#ifndef NDEBUG
//                    successors_[successors_.size() - 1] = nullptr;
//                    found = true;
//#endif
//                    successors_.pop_back();
//
//                    break;
//                }
//            }
//            assert(found);
//        }
//
//        void push_successor(Label* label)
//        {
//            successors_.push_back(label);
//        }
//
//        const std::vector<Label*>& successors()
//        {
//            return successors_;
//        }
//
//    private:
//        std::vector<Label*> successors_;
    };

    struct NodeEntry
    {
        using value_type = Label;

        bool in_queue;

        inline NodeEntry()
        :   in_queue(false),
            labels_(),
            labels_end_(0),
            new_labels_(),
            new_labels_end_(0)
        {
        }

        inline void push_back(Label* label)
        {
            add_label(new_labels_, new_labels_end_, label);
        }

        inline void erase(std::vector<Label*>& arr,
                   unsigned index)
        {

            assert(&arr == &labels_ ||
                   &arr == &new_labels_);

            delete arr[index];

            if(&arr == &labels_)
            {
                arr[index] = arr[labels_end_ - 1];
#ifndef NDEBUG
                arr[labels_end_ - 1] = nullptr;
#endif
                --labels_end_;
            }
            else
            {
                arr[index] = arr[new_labels_end_ - 1];
#ifndef NDEBUG
                arr[new_labels_end_ - 1] = nullptr;
#endif
                --new_labels_end_;
                
            }
        }

        inline void proceed_labels_it()
        {
            for(unsigned i = 0; i < new_labels_end_; ++i)
            {
                add_label(labels_, labels_end_, new_labels_[i]);
#ifndef NDEBUG
                new_labels_[i] = nullptr;
#endif
            }
            new_labels_end_ = 0;
        }

        inline bool has_new_labels()
        {
            return new_labels_end_ > 0;
        }

        inline std::vector<Label*>& labels()
        {
            return labels_;
        }

        inline std::vector<Label*>& new_labels()
        {
            return new_labels_;
        }

        inline unsigned labels_end()
        {
            return labels_end_;
        }

        inline unsigned new_labels_end()
        {
            return new_labels_end_;
        }

        ~NodeEntry()
        {
            for(unsigned i = 0; i < labels_end_; ++i)
            {
                delete labels_[i];
                labels_[i] = nullptr;
            }

            for(unsigned i = 0; i < new_labels_end_; ++i)
            {
                delete new_labels_[i];
                new_labels_[i] = nullptr;
            }
        }

    private:
        std::vector<Label*> labels_;
        unsigned labels_end_;

        std::vector<Label*> new_labels_;
        unsigned new_labels_end_;

        inline void add_label(std::vector<Label*>& arr,
                       unsigned& end_pointer,
                       Label* label)
        {

            if(arr.size() > end_pointer)
            {
                arr[end_pointer] = label;
            }
            else
            {
                arr.push_back(label);
            }
            ++end_pointer;
            
        }
    };

    bool check_domination(std::vector<Label*>& new_labels,
                          unsigned new_labels_end,
                          NodeEntry& neighbor_entry);

    void recursive_delete(Label& label);

    bool check_heuristic_prunable(const Label& label);


};

template<typename T>
class ring_buffer : private std::vector<T>
{
public:
    inline ring_buffer(unsigned capacity)
    :   std::vector<T>(capacity)
    {
    }

    inline void push_back(T object)
    {
        std::vector<T>::operator[](back_) = object;

        back_ = (back_ + 1) % std::vector<T>::size();
    }

    inline T front()
    {
        return std::vector<T>::operator[](front_);
    }

    inline void pop_front()
    {

#ifndef NDEBUG
        std::vector<T>::operator[](front_) = nullptr;
#endif

        front_ = (front_ + 1) % std::vector<T>::size();
    }

    inline bool empty()
    {
        return front_ == back_;
    }

private:
    unsigned front_ = 0;
    unsigned back_ = 0;
};

}

#endif /* BSSSA_H_ */
