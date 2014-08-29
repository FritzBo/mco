//
//  instance_scalarizer.h
//  mco
//
//  Created by Fritz Bökler on 29.08.14.
//
//

#ifndef mco_instance_scalarizer_h
#define mco_instance_scalarizer_h

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

class InstanceScalarizer {
public:
    static void scaleround_instance(const ogdf::Graph& graph,
                                    const ogdf::EdgeArray<Point>& cost_array,
                                    unsigned dimension,
                                    const Point& factor,
                                    ogdf::EdgeArray<Point>& new_cost_array);
};

}

#endif
