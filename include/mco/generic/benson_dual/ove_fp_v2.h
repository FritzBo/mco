//
//  ove_fp_v2.h
//  mco
//
//  Created by Fritz Bökler on 08.04.14.
//
//

#ifndef __mco__ove_fp_v2__
#define __mco__ove_fp_v2__

#include <queue>
#include <vector>
#include <algorithm>
#include <set>

#include <mco/generic/benson_dual/abstract_online_vertex_enumerator.h>
#include <mco/geometric/projective_geometry_utilities.h>

namespace mco {
    
class GraphlessOVE :
public AbstractOnlineVertexEnumerator {
public:
    virtual ~GraphlessOVE() {}
    
    GraphlessOVE(const Point& initial_value,
                 unsigned dimension,
                 double epsilon);
    
    template<typename ConstIterator>
    GraphlessOVE(unsigned dimension,
                 ConstIterator extreme_points_begin,
                 ConstIterator extreme_points_end,
                 ConstIterator extreme_rays_begin,
                 ConstIterator extreme_rays_end,
                 ConstIterator inequalities_begin,
                 ConstIterator inequalities_end,
                 double epsilon = 1E-9);
    
    inline bool has_next();
    
    inline Point * next_vertex();
    
    void add_hyperplane(Point &vertex, Point &normal, double rhs);
    
    unsigned int number_of_hyperplanes() { return inequalities_.size(); }

private:
    class GraphlessPoint : public Point {
    public:
        GraphlessPoint(unsigned dimension)
        :   Point(dimension) {}
        
        GraphlessPoint(Point&& that)
        :   Point(std::move(that)) {}
        
        std::list<unsigned> active_inequalities_;
        unsigned birth_index_;
        double* father_point_ = nullptr;
        bool removed = false;
    };
    
    bool check_adjacent(GraphlessPoint& point1,
                        const GraphlessPoint& point2);
    
    GraphlessPoint add_cut_point(const GraphlessPoint& outside_point,
                                 GraphlessPoint& inside_point,
                                 const Point& inequality);
    
    std::vector<GraphlessPoint> pending_points_;
    inline void push_pending(GraphlessPoint&& point);
    inline void pop_pending();
    inline GraphlessPoint& top_pending();
    
    std::list<GraphlessPoint> candidate_points_;
    
    std::list<GraphlessPoint> permanent_points_;
    std::vector<Point> inequalities_;
};
    
template<typename ConstIterator>
GraphlessOVE::GraphlessOVE(unsigned dimension,
                           ConstIterator extreme_points_begin,
                           ConstIterator extreme_points_end,
                           ConstIterator extreme_rays_begin,
                           ConstIterator extreme_rays_end,
                           ConstIterator inequalities_begin,
                           ConstIterator inequalities_end,
                           double epsilon)
: AbstractOnlineVertexEnumerator(dimension, epsilon) {
    
    for(auto it = extreme_points_begin; it != extreme_points_end; ++it) {
        Point new_point(dimension_ + 1);
        std::copy(it->cbegin(), it->cend(), new_point.begin());
        new_point[dimension_] = 1;
        pending_points_.push_back(std::move(new_point));
    }
    
    make_heap(pending_points_.begin(),
              pending_points_.end(),
              LexPointComparator());
    
    for(auto it = extreme_points_begin; it != extreme_points_end; ++it) {
        Point new_ray(dimension_ + 1);
        std::copy(it->cbegin(), it->cend(), new_ray.begin());
        new_ray[dimension_] = 0;
        pending_points_.push_back(std::move(new_ray));
    }
    
    std::copy(inequalities_begin,
              inequalities_end,
              back_inserter(inequalities_));
}
    
inline Point * GraphlessOVE::
next_vertex() {
    
    assert(is_heap(pending_points_.begin(),
                   pending_points_.begin() + pending_points_.size(),
                   LexPointComparator()));
    
    while(top_pending().removed) {
        pop_pending();
    }
    
#ifndef NDEBUG
    const double* adress = pending_points_.front().cbegin();
#endif
    
    candidate_points_.push_back(std::move(top_pending()));
    
    assert(candidate_points_.back().begin() == adress);
    
    pop_pending();
    
    return new Point(std::move(ProjectiveGeometry::from_projective<Point>(candidate_points_.back())));
}
    
inline bool GraphlessOVE::
has_next() {
    assert(is_heap(pending_points_.begin(),
                   pending_points_.end(),
                   LexPointComparator(epsilon_)));

    while(!pending_points_.empty() && top_pending().removed) {
        pop_pending();
    }
    
    return !pending_points_.empty();
}
    
inline void GraphlessOVE::
push_pending(GraphlessPoint&& point) {
    
#ifndef NDEBUG
    const double* adress = point.cbegin();
#endif
    
    pending_points_.push_back(std::move(point));
    
    assert(pending_points_.back().cbegin() == adress);
    
    std::push_heap(pending_points_.begin(),
                   pending_points_.end(),
                   LexPointComparator(epsilon_));
}

inline void GraphlessOVE::
pop_pending() {
    
    std::pop_heap(pending_points_.begin(),
                  pending_points_.end(),
                  LexPointComparator(epsilon_));
    pending_points_.pop_back();
}

inline auto GraphlessOVE::
top_pending() -> GraphlessPoint& {
    return pending_points_.front();
}
    
}

#endif /* defined(__mco__ove_fp_v2__) */
