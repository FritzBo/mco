//
//  ep_ls_module.h
//  mco
//
//  Created by Fritz Bökler on 18.09.2016.
//
//

#ifndef __mco__ep_ls_module__
#define __mco__ep_ls_module__

#include <mco/basic/forward_star.h>

#include "../basic/modules.h"

class EpLsModule : public AlgorithmModule<std::list<mco::node>> {
    
public:
    EpLsModule()
    :   AlgorithmModule(FORWARDSTAR) { }

    virtual void perform(int argc, char** args);
    virtual ~EpLsModule() {}
    
    virtual const std::list<std::pair<const std::list<mco::node>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:

    // stats:
    unsigned        num_nodes_                          = 0;
    unsigned        num_edges_                          = 0;
    unsigned        num_objectives_                     = 0;
    unsigned long   label_compares_                     = 0;
    unsigned        arc_pushes_                         = 0;
    unsigned long   deleted_labels_                     = 0;

    double          solution_time_                      = 0;
    
    std::list<std::pair<const std::list<mco::node>, const mco::Point>> solutions_;
};

#endif /* defined(__mco__ep_bs_module__) */
