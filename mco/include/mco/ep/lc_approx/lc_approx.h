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

    void set_bound(Point bound) {
        bound_ = std::move(bound);
        use_bounds_ = true;
    }
    
    void add_disjunctive_bound(Point bounds) {
        disj_bounds_.push_back(std::move(bounds));
        use_bounds_ = true;
    }

    template<typename InputIterator>
    void add_disjunctive_bounds(InputIterator begin,
                                InputIterator end) {

        disj_bounds_.insert(disj_bounds_.end(),
                            begin, end);
        use_bounds_ = true;
    }

private:
    
    Point min_e_;
    Point epsilon_;
    unsigned dimension_;
    
    bool use_heuristic_ = false;
    heuristic_type heuristic_;
    bool use_bounds_ = false;
    std::list<Point> disj_bounds_;
    Point bound_;
    
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
            sum(0),
            pred_edge(p_edge),
            pred_label(p_label),
            deleted(false) {


            auto it = cost.cbegin();
            while(it != cost.cend()) {
                sum += *it;
                ++it;
            }


        }

        Label(Label&& other)
        :   cost(std::move(other.cost)),
            pos(std::move(other.pos)),
            sum(other.sum),
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
            labels_end_(0),
            new_labels_(),
            new_labels_end_(0)
            { }

        NodeEntry(const NodeEntry& other)
        :   in_queue(false),
            labels_(),
            labels_end_(0),
            new_labels_(),
            new_labels_end_(0)
            { }

        NodeEntry(NodeEntry&& other)
        :   in_queue(false),
            labels_(),
            labels_end_(0),
            new_labels_(),
            new_labels_end_(0)
        { }

        void push_back(Label* label) {
            add_label(new_labels_, new_labels_end_, label);
        }
        
        void erase(std::vector<Label*>& arr,
                   unsigned index) {

            assert(&arr == &labels_ ||
                   &arr == &new_labels_);

            delete arr[index];

            if(&arr == &labels_) {
                arr[index] = arr[labels_end_ - 1];
#ifndef NDEBUG
                arr[labels_end_ - 1] = nullptr;
#endif
                --labels_end_;
            } else {
                arr[index] = arr[new_labels_end_ - 1];
#ifndef NDEBUG
                arr[new_labels_end_ - 1] = nullptr;
#endif
                --new_labels_end_;

            }
        }
        
        void proceed_labels_it() {
            for(unsigned i = 0; i < new_labels_end_; ++i) {
                add_label(labels_, labels_end_, new_labels_[i]);
#ifndef NDEBUG
                new_labels_[i] = nullptr;
#endif
            }
            new_labels_end_ = 0;
        }
        
        bool has_new_labels() {
            return new_labels_end_ > 0;
        }
        
        std::vector<Label*>& labels() {
            return labels_;
        }

        std::vector<Label*>& new_labels() {
            return new_labels_;
        }

        unsigned labels_end() {
            return labels_end_;
        }

        unsigned new_labels_end() {
            return new_labels_end_;
        }

    private:
        std::vector<Label*> labels_;
        unsigned labels_end_;

        std::vector<Label*> new_labels_;
        unsigned new_labels_end_;

        void add_label(std::vector<Label*>& arr,
                       unsigned& end_pointer,
                       Label* label) {

//            cout << arr.size() << endl;
            if(arr.size() > end_pointer) {
//                cout << "set" << endl;
                arr[end_pointer] = label;
            } else {
                arr.push_back(label);
//                cout << "push" << endl;
            }
            ++end_pointer;

        }
    };
    
    bool check_domination(std::vector<Label*>& new_labels,
                          unsigned new_labels_end,
                          NodeEntry& neighbor_entry);
    
    bool check_heuristic_prunable(const Label& label,
                                  const Point bound,
                                  const std::list<Point>& bounds);

    void recursive_delete(Label& label);
    
};
    
} // namespace mco

#endif /* defined(__mco__lc_approx__) */
