//
//  lc_approx.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__lc_approx__
#define __mco__lc_approx__

#include <ogdf/basic/Graph.H>

#include <mco/basic/abstract_solver.h>

namespace mco {

class LCApprox : public AbstractSolver<std::list<ogdf::edge>> {
    using cost_function_type = std::function<Point&(ogdf::edge)>;
public:
    
    void Solve(const ogdf::Graph& graph,
               cost_function_type cost_function,
               unsigned dimension,
               const ogdf::node source,
               const ogdf::node target,
               bool directed,
               double epsilon);
    
private:
    
    Point min_e_;
    double epsilon_;
    unsigned dimension_;
    
    void compute_pos(const Point& cost, Point& pos) const {
        for(unsigned i = 0; i < dimension_; ++i) {
            if(min_e_[i] != 0 && cost[i] != 0) {
                pos[i] = floor(log(cost[i] / min_e_[i]) / log(epsilon_));
            } else if(cost[i] != 0){
                pos[i] = floor(log(cost[i]) / log(epsilon_));
            } else {
                pos[i] = 0;
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
        const ogdf::node n;

        Label(const Point cost,
              const ogdf::node n,
              const LCApprox& app)
        :   cost(cost),
            pos(app.compute_pos(cost)),
            n(n) { }
    };
    
    struct NodeEntry {
        using value_type = Label;
        
        NodeEntry()
        :   labels_(),
            labels_it_(labels_.begin()) { }
        
        void push_back(Label&& label) {
            if(labels_it_ == labels_.end()) {
                labels_it_ = labels_.insert(labels_it_, std::move(label));
            } else {
                labels_.push_back(std::move(label));
            }
        }
        
        const std::list<Label>::const_iterator erase(std::list<Label>::const_iterator it) {
            if(it != labels_it_) {
                return labels_.erase(it);
            } else {
                return labels_it_ = labels_.erase(it);
            }
        }
        
        const std::list<Label>::const_iterator labels_it() {
            return labels_it_;
        }
        
        void proceed_labels_it() {
            labels_it_ = labels_.end();
        }
        
        bool has_new_labels() {
            return labels_it_ != labels_.end();
        }
        
        const std::list<Label>& labels() {
            return labels_;
        }
    private:
        std::list<Label> labels_;
        std::list<Label>::const_iterator labels_it_;
    };
    
    bool check_domination(std::list<Label>& new_labels,
                          NodeEntry& neighbor_entry);
    
};
    
} // namespace mco

#endif /* defined(__mco__lc_approx__) */
