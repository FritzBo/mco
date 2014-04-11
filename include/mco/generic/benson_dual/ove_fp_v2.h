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

#include <mco/generic/benson_dual/abstract_online_vertex_enumerator.h>

namespace mco {

class GraphLessOVE :
public AbstractOnlineVertexEnumerator {
public:
    virtual ~GraphLessOVE();
    
    GraphLessOVE(const Point& initial_value, unsigned dimension, double epsilon);
    
    bool has_next() { return !pending_points_.empty(); }
    
    Point * next_vertex();
    
    void add_hyperplane(Point &vertex, Point &normal, double rhs);
    
    unsigned int number_of_hyperplanes() { return hyperplanes_.size(); }

private:
    std::priority_queue<Point, std::vector<Point>, LexPointComparator> pending_points_;
    std::list<Point> permanent_points_;
    
    std::list<Point> hyperplanes_;
};
    
}

#endif /* defined(__mco__ove_fp_v2__) */
