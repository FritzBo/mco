//
//  ep_lc_approx_module.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__ep_lc_approx_module__
#define __mco__ep_lc_approx_module__

#include <mco/basic/forward_star.h>

#include "../basic/modules.h"

class EpLCApproxModule : public AlgorithmModule<std::list<mco::node>> {
    
public:
    EpLCApproxModule()
    :   AlgorithmModule(FORWARDSTAR) { }

    virtual void perform(int argc, char** args);
    virtual ~EpLCApproxModule() {}
    
    virtual const std::list<std::pair<const std::list<mco::node>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<mco::node>, const mco::Point>> solutions_;

    void parse_epsilon(const TCLAP::MultiArg<std::string>& epsilon_argument,
                       unsigned dimension,
                       mco::Point& epsilon);
    
};

#endif /* defined(__mco__ep_lc_approx_module__) */
