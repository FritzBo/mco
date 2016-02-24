//
//  paretoprep.h
//  mco
//
//  Created by Fritz Bökler on 18.02.16.
//
//

#ifndef __mco__paretoprep__
#define __mco__paretoprep__

#include <mco/basic/point.h>
#include <mco/basic/forward_star.h>

namespace mco {

class ParetoPrep {

public:
    void preprocess(const ForwardStar& reverse_star,
                    const ReverseEdgeArray<Point>& weights,
                    unsigned no_objectives,
                    node source,
                    node target);

private:
    class NodeEntry
    {
    public:
        Point temp_lower_bounds;
        std::vector<edge> predecessor_edges;

        NodeEntry(unsigned no_objectives)
        :   temp_lower_bounds(no_objectives),
            predecessor_edges(no_objectives)
        { }

    };

    void construct_path(const ForwardStar& reverse_star,
                        const ReverseEdgeArray<Point>& weights,
                        const FSNodeArray<NodeEntry>& node_entries,
                        node source,
                        node target,
                        unsigned objective,
                        Point& new_upper_bound);
};
}

#endif /* defined(__mco__paretoprep__) */
