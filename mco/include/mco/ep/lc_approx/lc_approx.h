//
//  lc_approx.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__lc_approx__
#define __mco__lc_approx__

#include <mco/basic/forward_star.h>
#include <mco/basic/abstract_solver.h>

namespace mco {

class LCApprox : public AbstractSolver<std::list<node>> {
    using cost_function_type = std::function<Point&(edge)>;
public:

    struct InstanceDescription
    {
        ForwardStar* graph;
        cost_function_type cost_function;
        unsigned dimension;
        node source;
        node target;
        Point* epsilon;
    };

    void Solve(const InstanceDescription& inst_desc);

private:
    
    Point min_e_;
    Point epsilon_;
    unsigned dimension_;

    void compute_pos(const Point& cost, Point& pos) const {
        for(unsigned i = 0; i < dimension_; ++i) {
            if(epsilon_[i] != 0.0) {
                if(min_e_[i] != 0 && cost[i] != 0) {
                    pos[i] = ceil(log(cost[i] / min_e_[i]) / log(epsilon_[i]));
                } else if(cost[i] != 0){
                    pos[i] = ceil(log(cost[i]) / log(epsilon_[i]));
                } else {
                    pos[i] = 0;
                }
            } else {
                pos[i] = cost[i];
            }
        }
    }
    
    Point compute_pos(const Point& cost) const {
        Point pos(dimension_);
        compute_pos(cost, pos);
        return pos;
    }
    
    struct Label {
        const Point cost;
        const Point pos;
        double sum;
        const node n;
        const edge pred_edge;
        Label* pred_label;
        std::list<Label*> succ_label;
        bool deleted = false;

        Label(const Point cost,
              const node n,
              const edge p_edge,
              Label* p_label,
              const LCApprox& app)
        :   cost(cost),
            pos(app.compute_pos(cost)),
            sum(0),
            n(n),
            pred_edge(p_edge),
            pred_label(p_label),
            deleted(false) {

            auto it = cost.cbegin();
            while(it != cost.cend()) {
                sum += *it;
                ++it;
            }
        }
    };
    
    struct NodeEntry {
        using value_type = Label;
        
        bool in_queue;
        
        NodeEntry()
        :   in_queue(false),
            labels_(0),
            new_labels_(0)
        { }

        void push_back(Label* label)
        {
            new_labels_.push_back(label);
        }
        
        void erase(std::vector<Label*>& arr,
                   unsigned index) {

            assert(&arr == &labels_ ||
                   &arr == &new_labels_);

            delete arr[index];

            arr[index] = arr[arr.size() - 1];
#ifndef NDEBUG
            arr[arr.size() - 1] = nullptr;
#endif
            arr.pop_back();
        }
        
        void proceed_labels_it() {
            for(unsigned i = 0; i < new_labels_.size(); ++i) {
                labels_.push_back(new_labels_[i]);
#ifndef NDEBUG
                new_labels_[i] = nullptr;
#endif
            }

            new_labels_.clear();
        }
        
        bool has_new_labels() {
            return !new_labels_.empty();
        }
        
        std::vector<Label*>& labels() {
            return labels_;
        }

        std::vector<Label*>& new_labels() {
            return new_labels_;
        }

    private:
        std::vector<Label*> labels_;
        std::vector<Label*> new_labels_;
    };
    
    bool check_domination(std::vector<Label*>& new_labels,
                          NodeEntry& neighbor_entry);
    
    void recursive_delete(Label& label);
    
};
    
} // namespace mco

#endif /* defined(__mco__lc_approx__) */
