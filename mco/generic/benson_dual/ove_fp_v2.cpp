//
//  ove_fp_v2.cpp
//  mco
//
//  Created by Fritz Bökler on 08.04.14.
//
//

#include <mco/generic/benson_dual/ove_fp_v2.h>

#include <list>

using std::list;

#include <mco/basic/point.h>

namespace mco {

GraphlessOVE::
GraphlessOVE(const Point& initial_value, unsigned dimension, double epsilon)
:   AbstractOnlineVertexEnumerator(dimension, epsilon),
    pending_points_(LexPointComparator()) {
        
    for(unsigned int i = 0; i < dimension_ - 1; ++i) {
        GraphlessPoint new_extreme_point(dimension_ + 1);
        
        for(unsigned int j = 0; j < dimension_ - 1; ++j) {
            new_extreme_point[j] = i == j ? 1 : 0;
            
            if(i != j) {
                new_extreme_point.active_inequalities_.push_back(j);
            }
        }
        
        new_extreme_point[dimension_ - 1] = initial_value[i];
        new_extreme_point[dimension_] = 1;
        
        new_extreme_point.active_inequalities_.push_back(dimension_ - 1);
        new_extreme_point.active_inequalities_.push_back(dimension_);
        new_extreme_point.birth_index_ = dimension_;

        pending_points_.push(std::move(new_extreme_point));
        
        Point new_inequality(dimension_ + 1);
        for(unsigned int j = 0; j < dimension_; ++j)
            new_inequality[j] = i == j ? 1 : 0;
        new_inequality[dimension_] = 0;
        
        inequalities_.push_back(std::move(new_inequality));
        
    }
    
    Point new_inequality(dimension_ + 1);
        
    for(unsigned int j = 0; j < dimension_ - 1; ++j) {
        new_inequality[j] = -1;
    }
        
    new_inequality[dimension_ - 1] = 0;
    new_inequality[dimension_] = 1;
    
    inequalities_.push_back(std::move(new_inequality));
    
    GraphlessPoint new_extreme_point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_; ++j) {
        new_extreme_point[j] = 0;
    }
    
    new_extreme_point[dimension_ - 1] = initial_value[dimension_ - 1];
    new_extreme_point[dimension_] = 1;
    
    for(unsigned int i = 0; i < dimension_ - 1; ++i) {
        new_extreme_point.active_inequalities_.push_back(i);
    }
    
    new_extreme_point.active_inequalities_.push_back(dimension_);
    new_extreme_point.birth_index_ = dimension_;
    
    pending_points_.push(std::move(new_extreme_point));
    
    
    new_inequality = Point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_; ++j) {
        new_inequality[j] = initial_value[j] - initial_value[dimension_ - 1];
    }
    new_inequality[dimension_ - 1] = -1;
    new_inequality[dimension_] = initial_value[dimension_ - 1];
    
    inequalities_.push_back(std::move(new_inequality));
    
    new_extreme_point = GraphlessPoint(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_ - 1; ++j) {
        new_extreme_point[j] = 0;
    }
    new_extreme_point[dimension_ - 1] = -1;
    new_extreme_point[dimension_] = 0;
    
    for(unsigned int i = 0; i < dimension_; ++i) {
        new_extreme_point.active_inequalities_.push_back(i);
    }
    
    new_extreme_point.birth_index_ = dimension_ - 1;
}
    
void GraphlessOVE::
add_hyperplane(mco::Point &vertex, mco::Point &normal, double rhs) {
    
}

}