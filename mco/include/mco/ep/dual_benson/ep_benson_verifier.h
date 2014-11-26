//
//  ep_benson_verifier.h
//  mco
//
//  Created by Fritz Bökler on 11.11.14.
//
//

#ifndef mco_ep_benson_verifier_h
#define mco_ep_benson_verifier_h

#include <mco/generic/benson_dual/upper_image_container.h>

class EpUpperImageVerifier {
public:
    template<typename OutputIterator>
    void verify(AbstractUpperImageContainer<OutputIterator>& container);

private:
    void verify_double_description();
    void verify_extreme_points();
    void verify_inequalities();
};

template<typename OutputIterator>
void EpUpperImageVerifier::
verify(AbstractUpperImageContainer<OutputIterator>& container) {

}

#endif
