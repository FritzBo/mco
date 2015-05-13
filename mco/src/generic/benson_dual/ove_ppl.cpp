/*
 * online_vertex_enumerator.cpp
 *
 *  Created on: 30.09.2013
 *      Author: fritz
 */

#include <mco/generic/benson_dual/ove_ppl.h>

#include <iostream>
#include <cassert>
#include <list>

using std::list;

using PPL::Constraint_System;
using PPL::Linear_Expression;
using PPL::NNC_Polyhedron;
using PPL::Generator_System;
using PPL::Generator;
using PPL::Coefficient;

namespace mco {

OnlineVertexEnumeratorPPL::
OnlineVertexEnumeratorPPL(Point &initial_value,
                          unsigned int dimension,
                          double epsilon)
    :   AbstractOnlineVertexEnumerator(dimension, epsilon),
        number_hyperplanes_(0)

{

	for(unsigned int i = 0; i < dimension_ - 1; ++i) {
        h_representation_.insert(PPL::Variable(i) >= 0);

		Point p(dimension_);
		for(unsigned int j = 0; j < dimension_ - 1; ++j)
			p[j] = i == j ? 1 : 0;
		p[dimension_ - 1] = initial_value[i];

		unprocessed_vertices_.push_back(p);
	}

	Point p(dimension);
	p[dimension_ - 1] = initial_value[dimension_ - 1];

	unprocessed_vertices_.push_back(p);

    Linear_Expression lambda_upper_bound;
	for(unsigned int j = 0; j < dimension_ - 1; ++j) {
        lambda_upper_bound = PPL::Variable(j);
    }

    h_representation_.insert(lambda_upper_bound <= 1);

    Linear_Expression initial_inequality;

	for(unsigned int j = 0; j < dimension_ - 1; ++j) {
        initial_inequality += PPL::Variable(j) * (initial_value[j] - initial_value[dimension - 1]);
    }

    initial_inequality += PPL::Variable(dimension - 1) * -1;

	h_representation_.insert(initial_inequality >= - initial_value[dimension_ - 1]);

}

OnlineVertexEnumeratorPPL::~OnlineVertexEnumeratorPPL() {

}

bool OnlineVertexEnumeratorPPL::has_next() {
	return !unprocessed_vertices_.empty();
}

Point * OnlineVertexEnumeratorPPL::next_vertex() {
	if(unprocessed_vertices_.empty())
		return nullptr;

	Point * p = new Point(unprocessed_vertices_.front());
	unprocessed_vertices_.pop_front();

	return p;
}

void OnlineVertexEnumeratorPPL::add_hyperplane(Point &vertex, Point &normal, double rhs) {
	clock_t start = clock();
//	std::cout << start << std::endl;

//	std::cout << "begin" << std::endl;

//	std::cout << "normal: " << normal << std::endl;
//	std::cout << "rhs: " << rhs << std::endl;

	number_hyperplanes_++;

//    if(new_line_ == dimension_ + 1) {
//        cout << ddf_almostzero << endl;
//    }

    Linear_Expression new_inequality;
	for(unsigned int i = 0; i < dimension_; ++i) {
        new_inequality += normal[i] * PPL::Variable(i);
    }

    Constraint_System representation_copy(h_representation_);

    representation_copy.insert(new_inequality == rhs);
    h_representation_.insert(new_inequality >= rhs);

    NNC_Polyhedron poly(representation_copy, PPL::Recycle_Input());

	auto& generators = poly.generators();

	auto it = unprocessed_vertices_.begin();
	while(it != unprocessed_vertices_.end()) {
		if(*it * normal - rhs < epsilon_) {
//			std::cout << "Removed " << *it << std::endl;
			it = unprocessed_vertices_.erase(it);
		} else
			it++;
	}

//	std::cout << "adding.." << std::endl;

    for (Generator_System::const_iterator i = generators.begin(), gs_end = generators.end(); i != gs_end; ++i) {

        const Generator& g = *i;

        if(!g.is_point()) {
			continue;
        }

        const Coefficient& divisor = g.divisor();

		Point p(dimension_);
        for (unsigned j = 0; j < dimension_; ++j) {
            mpz_class numer, denom;

            PPL::assign_r(numer,
                          g.coefficient(PPL::Variable(j)),
                          PPL::ROUND_NOT_NEEDED);

            PPL::assign_r(denom,
                          divisor,
                          PPL::ROUND_NOT_NEEDED);

            p[j] = mpq_class(numer, denom).get_d();
        }

		unprocessed_vertices_.push_back(p);

//		std::cout << "Added " << p << std::endl;
	}

	cycles_ += clock() - start;
}

} /* namespace mco */
