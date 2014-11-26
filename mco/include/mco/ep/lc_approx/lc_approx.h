//
//  lc_approx.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__lc_approx__
#define __mco__lc_approx__

#include <ogdf/basic/Graph.h>

#include <mco/basic/abstract_solver.h>

namespace mco {

class LCApprox : public AbstractSolver<std::list<ogdf::edge>> {
    using cost_function_type = std::function<Point&(ogdf::edge)>;
    using heuristic_type = std::function<double(ogdf::node, unsigned)>;
public:
    
    void Solve(const ogdf::Graph& graph,
               cost_function_type cost_function,
               unsigned dimension,
               const ogdf::node source,
               const ogdf::node target,
               bool directed,
               const Point& epsilon);
    
    void set_heuristic(heuristic_type heuristic) {
        heuristic_ = heuristic;
        use_heuristic_ = true;
    }
    
    void set_bounds(Point bounds) {
        bounds_ = std::move(bounds);
        use_bounds_ = true;
    }
    
private:
    
    Point min_e_;
    Point epsilon_;
    unsigned dimension_;
    
    bool use_heuristic_ = false;
    heuristic_type heuristic_;
    bool use_bounds_ = false;
    Point bounds_;
    
    void compute_pos(const Point& cost, Point& pos) const {
        for(unsigned i = 0; i < dimension_; ++i) {
            if(epsilon_[i] != 0.0) {
                if(min_e_[i] != 0 && cost[i] != 0) {
                    pos[i] = floor(log(cost[i] / min_e_[i]) / log(epsilon_[i]));
                } else if(cost[i] != 0){
                    pos[i] = floor(log(cost[i]) / log(epsilon_[i]));
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
        const ogdf::node n;
        const ogdf::edge pred_edge;
        Label* pred_label;
        std::list<Label*> succ_label;
        bool deleted = false;

        Label(const Point cost,
              const ogdf::node n,
              const ogdf::edge p_edge,
              Label* p_label,
              const LCApprox& app)
        :   cost(cost),
            pos(app.compute_pos(cost)),
            n(n),
            pred_edge(p_edge),
            pred_label(p_label),
            deleted(false) { }

        Label(Label&& other)
        :   cost(std::move(other.cost)),
            pos(std::move(other.pos)),
            n(other.n),
            pred_edge(other.pred_edge),
            pred_label(other.pred_label) {
        }
    };
    
    struct NodeEntry {
        using value_type = Label;
        
        bool in_queue;
        
        NodeEntry()
        :   in_queue(false),
            labels_(),
            labels_it_(labels_.begin()) { }

        NodeEntry(const NodeEntry& other)
        :   in_queue(false),
            labels_(),
            labels_it_(labels_.begin()) {
        }

        NodeEntry(NodeEntry&& other)
        :   in_queue(false),
            labels_(),
            labels_it_(labels_.begin()) {

        }
        
        void push_back(Label&& label) {
            if(labels_it_ == labels_.end()) {
                labels_it_ = labels_.insert(labels_it_, std::move(label));
            } else {
                labels_.push_back(std::move(label));
            }
        }
        
        std::list<Label>::iterator erase(std::list<Label>::iterator it) {
            if(it != labels_it_) {
                return labels_.erase(it);
            } else {
                return labels_it_ = labels_.erase(it);
            }
        }
        
        const std::list<Label>::iterator labels_it() {
            return labels_it_;
        }
        
        void proceed_labels_it() {
            labels_it_ = labels_.end();
        }
        
        bool has_new_labels() {
            return labels_it_ != labels_.end();
        }
        
        std::list<Label>& labels() {
            return labels_;
        }
    private:
        std::list<Label> labels_;
        std::list<Label>::iterator labels_it_;
    };
    
    bool check_domination(std::list<Label>& new_labels,
                          NodeEntry& neighbor_entry);
    
    bool check_heuristic_prunable(const Label& label,
                                  const Point& bounds);

    void recursive_delete(Label& label);
    
};
    
} // namespace mco

#endif /* defined(__mco__lc_approx__) */
