//
//  ove_fp_v2.cpp
//  mco
//
//  Created by Fritz Bökler on 08.04.14.
//
//

#include <mco/generic/benson_dual/ove_fp_v2.h>

#include <list>
#include <utility>
#include <set>
#include <cmath>

using std::list;
using std::pair;
using std::make_pair;
using std::set;
using std::abs;
using std::make_heap;

#include <mco/basic/point.h>

namespace mco {
    
GraphlessOVE::
GraphlessOVE(const Point& initial_value, unsigned dimension, double epsilon)
:   AbstractOnlineVertexEnumerator(dimension, epsilon) {
        
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

        push_pending(std::move(new_extreme_point));
        
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
    
    push_pending(std::move(new_extreme_point));
    
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
    
    permanent_points_.push_back(std::move(new_extreme_point));
}
    
void GraphlessOVE::
add_hyperplane(Point &vertex, Point &normal, double rhs) {
    
    // Create projective version of the given hyperplane
    Point projective_hyperplane(dimension_ + 1);
    std::copy(normal.cbegin(), normal.cend(), projective_hyperplane.begin());
    projective_hyperplane[normal.dimension()] = -rhs;
    
    // All candidates but the last are permanent points
    auto end_it = candidate_points_.end();
    move(candidate_points_.begin(), --end_it, back_inserter(permanent_points_));
    
    // three sets of extreme points:
    // cutted points, on-plane points and inside points
    list<GraphlessPoint*> inside_points;
    list<GraphlessPoint*> cut_off_points;
    
    cut_off_points.push_back(&candidate_points_.back());
    
    // Check for each pending point to which set it belongs
    for(auto& point : pending_points_) {
        double distance = projective_hyperplane * point;
        
        if(distance < -epsilon_) {
            cut_off_points.push_back(&point);
        } else if(distance > epsilon_) {
            inside_points.push_back(&point);
        }
    }
    
    list<GraphlessPoint> new_points;
    
    // Check for each cutted off pending point if there is a permanent
    // or another pending point strictly inside the polytope which makes them a
    // candidate pair
    for(auto pending_it = cut_off_points.begin();
        pending_it != cut_off_points.end();
        ++pending_it) {
        
        GraphlessPoint& pend_point = **pending_it;
        
        if(pend_point.removed) {
            continue;
        }
        
        for(auto& perm_point : permanent_points_) {
            if(abs(perm_point * projective_hyperplane) > epsilon_) {
                if(check_adjacent(pend_point, perm_point)) {
                    
                    new_points.push_back(add_cut_point(pend_point,
                                                       perm_point,
                                                       projective_hyperplane)
                                         );
                }
            }
        }
        
        for(auto pending_it2 = inside_points.begin();
            pending_it2 != inside_points.end();
            ++pending_it2) {
            
            const GraphlessPoint& pend_point2 = **pending_it2;
            
            if(!pend_point2.removed) {
                if(check_adjacent(pend_point, pend_point2)) {
                    
                    new_points.push_back(add_cut_point(pend_point,
                                                       pend_point2,
                                                       projective_hyperplane)
                                         );
                }
            }
        }
    }
    
    // mark all cut off points as removed
    for(auto point : cut_off_points) {
        point->removed = true;
    }
    
    // put all new points in the priority queue
    for(auto& point : new_points) {
        push_pending(std::move(point));
    }
    
    candidate_points_.clear();
}
    
bool GraphlessOVE::
check_adjacent(GraphlessPoint& p1, const GraphlessPoint& p2) {
    
    if(abs(p1[dimension_]) < epsilon_ && abs(p2[dimension_]) < epsilon_) {
        return false;
    }
    
    if(p1.father_point_ == &p2 || p2.father_point_ == &p1) {
        return true;
    }
    
	if(dimension_ == 2) {
		return true;
    
	} else if(dimension_ == 3) {
        
		list<int> inequality_intersection;
                
		set_intersection(p1.active_inequalities_.cbegin(),
                         p1.active_inequalities_.cend(),
                         p2.active_inequalities_.cbegin(),
                         p2.active_inequalities_.cend(),
                         back_inserter(inequality_intersection));
        
		unsigned tight_inequalities = inequality_intersection.size();
        
		if(tight_inequalities == dimension_ - 1)
			return true;
		else if(tight_inequalities < dimension_ - 1)
			return false;
		else
			assert(false);
        
	} else {
		unsigned int i = 0;
		list<unsigned int> tight_inequalities;
        
		set_intersection(p1.active_inequalities_.cbegin(),
                         p1.active_inequalities_.cend(),
                         p2.active_inequalities_.cbegin(),
                         p2.active_inequalities_.cend(),
                         back_inserter(tight_inequalities));
        
		unsigned int k = max(p1.birth_index_, p2.birth_index_);
		bool nc2 = false;
        
		for(auto index: tight_inequalities) {
			if(index >= k && index <= inequalities_.size() - 1) {
				nc2 = true;
				break;
			}
        }
        
		// [FP96] NC2
		if(!nc2)
			return false;
        
		// [FP96] NC1
		if(tight_inequalities.size() < dimension_ - 2)
			return false;
        
        set<GraphlessPoint*, LexPointComparator>
        common_vertices((LexPointComparator()));
        
		for(auto inequality_index : tight_inequalities) {
            
			Point& inequality = inequalities_[inequality_index];
            
			set<GraphlessPoint *, LexPointComparator>
            new_vertices((LexPointComparator()));
            
			for(auto& test_point : permanent_points_) {
//				cout << "checking node " << n << " with point " << *node_points_[n] << ": " << (*inequality) * (*node_points_[n]) << endl;
				if(abs(inequality * test_point) < epsilon_) {
					new_vertices.insert(&test_point);
                }
			}
            
			if(common_vertices.empty()) {
				common_vertices.insert(new_vertices.begin(),
                                       new_vertices.end());
			} else {
				list<GraphlessPoint *> temp_points;
				set_intersection(common_vertices.begin(),
                                 common_vertices.end(),
                                 new_vertices.begin(),
                                 new_vertices.end(),
                                 back_inserter(temp_points),
                                 LexPointComparator());
                
//				cout << "intersection size: "<< temp_points.size() << endl;
				common_vertices.clear();
				common_vertices.insert(temp_points.begin(), temp_points.end());
			}
            
//			cout << "current number of common vertices: " << common_vertices.size() << endl;
            
			i += 1;
		}
        
        //		cout << "Number of common vertices: " << common_vertices.size() << endl;
        
		if(common_vertices.size() == 2)
			return true;
		else
			return false;
        
	}
    
	return false;
}
    
auto GraphlessOVE::
add_cut_point(const GraphlessPoint& point1,
              const GraphlessPoint& point2,
              const Point& inequality) -> GraphlessPoint {
    
    double alpha = - inequality * point1;
    GraphlessPoint diff_direction = point2 - point1;
    
    assert(abs(inequality * diff_direction) > epsilon_);
    
    alpha /=  inequality * diff_direction;
    
    assert(alpha > epsilon_ && alpha < 1 + epsilon_);
    
    diff_direction *= alpha;
    GraphlessPoint cut_point = point1 + diff_direction;
    ProjectiveGeometry::normalize_projective(cut_point);
    
    return cut_point;
}

}