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
using std::make_heap;
using std::vector;

#include <mco/basic/point.h>

namespace mco {
    
GraphlessOVE::
GraphlessOVE(const Point& initial_value, unsigned dimension, double epsilon)
:   AbstractOnlineVertexEnumerator(dimension, epsilon) {
    
    pending_points_.clear();
    
    GraphlessPoint* new_extreme_point;
    for(unsigned int i = 0; i < dimension_ - 1; ++i) {
        new_extreme_point = new GraphlessPoint(dimension_ + 1);
        
        for(unsigned int j = 0; j < dimension_ - 1; ++j) {
            new_extreme_point->operator[](j) = i == j ? 1 : 0;
            
            if(i != j) {
                new_extreme_point->active_inequalities_.push_back(j);
            }
        }
        
        new_extreme_point->operator[](dimension_ - 1) = initial_value[i];
        new_extreme_point->operator[](dimension_) = 1;
        
        new_extreme_point->active_inequalities_.push_back(dimension_ - 1);
        new_extreme_point->active_inequalities_.push_back(dimension_);
        new_extreme_point->birth_index_ = dimension_;

        push_pending(new_extreme_point);
        
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
    
    new_extreme_point = new GraphlessPoint(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_; ++j) {
        new_extreme_point->operator[](j) = 0;
    }
    
    new_extreme_point->operator[](dimension_ - 1) = initial_value[dimension_ - 1];
    new_extreme_point->operator[](dimension_) = 1;
    
    for(unsigned int i = 0; i < dimension_ - 1; ++i) {
        new_extreme_point->active_inequalities_.push_back(i);
    }
    
    new_extreme_point->active_inequalities_.push_back(dimension_);
    new_extreme_point->birth_index_ = dimension_;
    
    push_pending(new_extreme_point);
    
    new_inequality = Point(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_; ++j) {
        new_inequality[j] = initial_value[j] - initial_value[dimension_ - 1];
    }
    new_inequality[dimension_ - 1] = -1;
    new_inequality[dimension_] = initial_value[dimension_ - 1];
    
    inequalities_.push_back(std::move(new_inequality));
    
    
    // Point at infinity
    new_extreme_point = new GraphlessPoint(dimension_ + 1);
    for(unsigned int j = 0; j < dimension_ - 1; ++j) {
        new_extreme_point->operator[](j) = 0;
    }
    new_extreme_point->operator[](dimension_ - 1) = -1;
    new_extreme_point->operator[](dimension_) = 0;
    
    for(unsigned int i = 0; i < dimension_; ++i) {
        new_extreme_point->active_inequalities_.push_back(i);
    }
    
    new_extreme_point->birth_index_ = dimension_ - 1;
    new_extreme_point->set_father(new_extreme_point);
    
    permanent_points_.push_back(new_extreme_point);
    
    make_heap(pending_points_.begin(),
              pending_points_.end(),
              LexPointComparator(epsilon_));
    
    GraphlessPoint& infinity_point = *permanent_points_.front();
    
    for(GraphlessPoint* point : pending_points_) {
        point->set_father(&infinity_point);
    }
}
    
void GraphlessOVE::
add_hyperplane(Point &vertex, Point &normal, double rhs) {
    
    is_heap(pending_points_.begin(), pending_points_.end(),
            LexPointComparator(epsilon_));
    
    
    // Create projective version of the given hyperplane
    Point projective_hyperplane(dimension_ + 1);
    std::copy(normal.cbegin(), normal.cend(), projective_hyperplane.begin());
    projective_hyperplane[normal.dimension()] = -rhs;
    
    // Add new inequality
    inequalities_.push_back(projective_hyperplane);
    
    // All candidates but the last are permanent points
    GraphlessPoint* last_candidate = candidate_points_.back();
    for(auto point : candidate_points_) {
        if(point != last_candidate) {
            permanent_points_.push_back(point);
        }
    }
    candidate_points_.clear();
    
    // three sets of extreme points:
    // cutted points, on-plane points and inside points
    list<GraphlessPoint*> inside_points;
    list<GraphlessPoint*> cut_off_points;
    
    pending_points_.push_back(last_candidate);
    
    // Check for each pending point to which set it belongs
    for(auto point : pending_points_) {
        double distance = projective_hyperplane * *point;
        
        // point is cut off by the inequality
        if(distance < -epsilon_) {
            cut_off_points.push_back(point);
        
        // point is inside the inequality induced halfspace
        } else if(distance > epsilon_) {
            inside_points.push_back(point);
            
        // point is on the inequality induced hyperplane
        } else {
            point->active_inequalities_.push_back(inequalities_.size() - 1);
        }
    }
    
    // list to store new points
    list<GraphlessPoint*> new_points;
    
    // check for each cutted off pending point if there is a permanent
    // or another pending point strictly inside the polytope which makes them a
    // candidate pair
    for(auto pending_it = cut_off_points.begin();
        pending_it != cut_off_points.end();
        ++pending_it) {
        
        GraphlessPoint& pend_point = **pending_it;
        
        assert(pend_point.active_inequalities_.size() >= dimension_);
        
        if(pend_point.removed) {
            continue;
        }
        
        for(auto perm_point : permanent_points_) {
            double distance = *perm_point * projective_hyperplane;
            assert(perm_point->active_inequalities_.size() >= dimension_);
            
            assert(distance > -epsilon_);
            
            if(distance > epsilon_) {
                if(check_adjacent(pend_point, *perm_point)) {
                    
                    new_points.push_back(add_cut_point(pend_point,
                                                       *perm_point,
                                                       projective_hyperplane)
                                         );
                }
                
            } else {
                perm_point->active_inequalities_.push_back(inequalities_.size() - 1);
            }
        }
        
        for(auto pending_it2 = inside_points.begin();
            pending_it2 != inside_points.end();
            ++pending_it2) {
            
            GraphlessPoint& pend_point2 = **pending_it2;
            assert(pend_point2.active_inequalities_.size() >= dimension_);
            
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
    
    delete pending_points_.back();
    pending_points_.pop_back();
    
    // put all new points in the priority queue
    for(auto point : new_points) {
        push_pending(point);
    }
}
    
    
bool GraphlessOVE::
check_adjacent(GraphlessPoint& p1, const GraphlessPoint& p2) {
    
#ifndef NDEBUG
    bool debug_output = false;
    if(debug_output) {
        cout << "Checking adjacency of " << p1 << " and " << p2 << endl;
        cout << "With father points ";
        if(p1.father() != nullptr) {
            cout << *p1.father();
        } else {
            cout << "NULL";
        }
        cout << " and ";
        if(p2.father() != nullptr) {
            cout << *p2.father();
        } else {
            cout << "NULL";
        }
        cout << "." << endl;
    }
#endif
    
    assert(!EqualityPointComparator()(p1, p2));
    
    if(std::abs(p1[dimension_]) < epsilon_ && std::abs(p2[dimension_]) < epsilon_) {
#ifndef NDEBUG
        if(debug_output) {
            cout << "Since both points are at infinity, they cannot be adjacent." << endl;
        }
#endif
        return false;
    }
    
    if(p1.father() == &p2 || p2.father() == &p1) {
#ifndef NDEBUG
        if(debug_output) {
            cout << "One is the father of the other, so they are adjacent." << endl;
        }
#endif
        return true;
    }
    
#ifndef NDEBUG
    // Check, if one not is really not the father of the other
    bool equal = true;
    if(p1.father() != nullptr) {
        for(unsigned i = 0; i < dimension_ + 1; ++i) {
            if(std::abs(p1.father()->operator[](i) - p2[i]) > epsilon_) {
                equal = false;
                break;
            }
        }
        assert(!equal);
    }
    if(p2.father() != nullptr) {
        for(unsigned i = 0; i < dimension_ + 1; ++i) {
            if(std::abs(p2.father()->operator[](i) - p1[i]) > epsilon_) {
                equal = false;
                break;
            }
        }
        assert(!equal);
    }
#endif
    
    if(dimension_ <= 3) {
        
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
		if(!nc2) {
#ifndef NDEBUG
            if(debug_output) {
                cout << "Not adjacent, because auf NC2" << endl;
            }
#endif
			return false;
        }
        
		// [FP96] NC1
		if(tight_inequalities.size() < dimension_ - 2) {
#ifndef NDEBUG
            if(debug_output) {
                cout << "Not adjacent, because of NC1" << endl;
            }
#endif
			return false;
        }
        
        set<GraphlessPoint*> common_vertices;
        
        assert(candidate_points_.empty());
        
		for(auto inequality_index : tight_inequalities) {
            
			Point& inequality = inequalities_[inequality_index];
            
			set<GraphlessPoint *> new_vertices;
            
#ifndef NDEBUG
            if(debug_output) {
                cout << "Checking permanent points..." << endl;
            }
#endif
            
			for(auto test_point : permanent_points_) {
                
#ifndef NDEBUG
                if(debug_output) {
                    cout << "checking point " << *test_point << ": " <<
                    inequality * *test_point;
                }
#endif
                
				if(std::abs(inequality * *test_point) < epsilon_) {
					new_vertices.insert(test_point);
#ifndef NDEBUG
                    if(debug_output) {
                        cout << " X" << endl;
                    }
#endif
                }
#ifndef NDEBUG
                else {
                    if(debug_output) {
                        cout << endl;
                    }
                }
#endif
			}
            
#ifndef NDEBUG
            if(debug_output) {
                cout << "Checking pending..." << endl;
            }
#endif
            
            for(auto test_point : pending_points_) {
                
                if(!test_point->removed) {
#ifndef NDEBUG
                    if(debug_output) {
                        cout << "checking point " << *test_point << ": " <<
                        inequality * *test_point;
                    }
#endif
                    if(std::abs(inequality * *test_point) < epsilon_) {
                        new_vertices.insert(test_point);
#ifndef NDEBUG
                        if(debug_output) {
                            cout << " X" << endl;
                        }
#endif
                    }
#ifndef NDEBUG
                    else {
                        if(debug_output) {
                            cout << endl;
                        }
                    }
                    
#endif
                }   // removed check
                
			}   // pending points loop
            
			if(common_vertices.empty()) {
				common_vertices.insert(new_vertices.begin(),
                                       new_vertices.end());
			} else {
				list<GraphlessPoint *> temp_points;
				set_intersection(common_vertices.begin(),
                                 common_vertices.end(),
                                 new_vertices.begin(),
                                 new_vertices.end(),
                                 back_inserter(temp_points));
#ifndef NDEBUG
                if(debug_output) {
                    cout << "intersection size: "<< temp_points.size() << endl;
                }
#endif
				common_vertices.clear();
				common_vertices.insert(temp_points.begin(), temp_points.end());
			}
            
#ifndef NDEBUG
            if(debug_output) {
                cout << "current number of common points: " << common_vertices.size() << endl;
            }
#endif
            
			i += 1;

		}
        
#ifndef NDEBUG
        if(debug_output) {
            cout << "Number of common points: " << common_vertices.size() << endl;
        }
#endif
        
		if(common_vertices.size() == 2)
			return true;
		else
			return false;
        
	}
    
	return true;
}
    
auto GraphlessOVE::
add_cut_point(const GraphlessPoint& outside_point,
              GraphlessPoint& inside_point,
              const Point& inequality) -> GraphlessPoint* {
    
    
#ifndef NDEBUG
    bool debug_output = false;
    if(debug_output)
        cout << "+-> Adding new point between outside point " << outside_point <<
        " and inside point " << inside_point << " : ";
#endif
    
    double alpha = - inequality * outside_point;
    Point diff_direction = inside_point - outside_point;
    
    assert(std::abs(inequality * diff_direction) > epsilon_);
    
    alpha /=  inequality * diff_direction;
    
    assert(alpha > epsilon_ && alpha < 1 + epsilon_);
    
    diff_direction *= alpha;
    GraphlessPoint* cut_point = new GraphlessPoint(outside_point + diff_direction);
    ProjectiveGeometry::normalize_projective(*cut_point);
    
//    set_intersection(outside_point.active_inequalities_.cbegin(),
//                     outside_point.active_inequalities_.cend(),
//                     inside_point.active_inequalities_.cbegin(),
//                     inside_point.active_inequalities_.cend(),
//                     back_inserter(cut_point->active_inequalities_));
    
//    cut_point->active_inequalities_.push_back(inequalities_.size() - 1);
        
    unsigned i = 0;
    for(auto inequality : inequalities_) {
        double distance = *cut_point * inequality;
#ifndef NDBUG
        if(cut_point->operator[](0) >= 0.0001298 && cut_point->operator[](0) <= 0.0001299) {
            cout << "Inequality: " << inequality << ", distance: " <<
            distance << " " << (std::abs(distance) < epsilon_ ? "X" : "") << endl;
        }
#endif
        if(std::abs(distance) < epsilon_) {
            cut_point->active_inequalities_.push_back(i);
        }
        ++i;
    }
    
    cut_point->set_father(&inside_point);
    cut_point->birth_index_ = inequalities_.size() - 1;
        
#ifndef NDEBUG
    if(debug_output) {
        cout << *cut_point << endl;
    }
#endif

#ifndef NDEBUG
    if(cut_point->active_inequalities_.size() < dimension_) {
        cout << "New point " << *cut_point << " has only " <<
        cut_point->active_inequalities_.size() << " active inequalities" << endl;
        assert(false);
    }
#endif
    
    return cut_point;
}

}