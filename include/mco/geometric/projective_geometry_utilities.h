//
//  projective_geometry_utilities.h
//  mco
//
//  Created by Fritz Bökler on 11.04.14.
//
//

#ifndef __mco__projective_geometry_utilities__
#define __mco__projective_geometry_utilities__

#include <type_traits>

#include <mco/basic/point.h>

namespace mco {

class ProjectiveGeometry {
public:
    template<typename PointType>
    inline PointType to_projective(const Point& point);
    
    inline void normalize_projective(Point& point);
};

template<typename PointType,
         typename std::enable_if<std::is_base_of<Point, PointType>>::type>
inline PointType NodeListVE::
to_projective(const Point &point) {
    
    unsigned dimension = point.dimension();
    PointType new_proj_point(dimension + 1);
    
	for(unsigned int i = 0; i < dimension; ++i)
		new_proj_point[i] = point[i];
    
	new_proj_point[dimension_] = 1;
    
	return new_proj_point;
}
    
inline void ProjectiveGeometry::
normalize_projective(Point& projective_point) {
    
    for(unsigned int i = 0; i < dimension_; ++i)
        projective_point[i] = projective_point[i] / projective_point[dimension_];
    
    normalized_projective_values[dimension_] = 1;
}
    
}


#endif /* defined(__mco__projective_geometry_utilities__) */
