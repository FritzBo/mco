//
//  ep_lc_approx_module.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__ep_lc_approx_module__
#define __mco__ep_lc_approx_module__

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EpLCApproxModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~EpLCApproxModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;
    
    
};

#endif /* defined(__mco__ep_lc_approx_module__) */
