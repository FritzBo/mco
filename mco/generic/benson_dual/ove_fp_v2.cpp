//
//  ove_fp_v2.cpp
//  mco
//
//  Created by Fritz Bökler on 08.04.14.
//
//

#include <mco/generic/benson_dual/ove_fp_v2.h>

#include <mco/basic/point.h>

namespace mco {

GraphLessOVE::
GraphLessOVE(const Point& initial_value, unsigned dimension, double epsilon)
:   AbstractOnlineVertexEnumerator(dimension, epsilon),
    pending_points_(LexPointComparator()) {
    
}

}