#pragma once
/*
 * label.h
 *
 *  Created on: 27.03.2013
 *      Author: fritz
 */

#ifndef LABEL_H_
#define LABEL_H_

#include <ogdf/basic/Graph.h>

#include <mco/basic/point.h>

namespace mco {

struct Label {
	const Point point;
	ogdf::node n;
	const Label * const pred;
	bool mark_dominated;

    inline Label(Point&& point,
                 ogdf::node n,
                 const Label *pred)
    :   point(std::move(point)),
        n(n),
        pred(pred),
        mark_dominated(false)
    {
    }

    inline Label(const Label &label)
    :   point(label.point),
        n(label.n),
        pred(label.pred),
        mark_dominated(label.mark_dominated)
    {
    }

	Label & operator=(const Label &label) = delete;
};

}

#endif /* LABEL_H_ */
