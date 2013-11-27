#pragma once
/*
 * dual_benson_vertex_container.h
 *
 *  Created on: 01.10.2013
 *      Author: fritz
 */

#ifndef DUAL_BENSON_VERTEX_CONTAINER_CDD_H_
#define DUAL_BENSON_VERTEX_CONTAINER_CDD_H_

#include <mco/generic/geometric/online_vertex_enumerator_cdd.h>

namespace mco {

class DualBensonVertexContainerCDD : public OnlineVertexEnumeratorCDD {
public:
	DualBensonVertexContainerCDD(Point &initial_value, unsigned int dimension, double epsilon);
};

} /* namespace mco */
#endif /* DUAL_BENSON_VERTEX_CONTAINER_H_ */
