//
//  benson_module.h
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#ifndef __mco__ap_benson_module__
#define __mco__ap_benson_module__

#include <list>

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class ApBensonModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:

    ApBensonModule()
    :   AlgorithmModule(OGDF) { }

    virtual void perform(int argc, char** args);
    virtual ~ApBensonModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;

    unsigned objectives_;
    unsigned nodes_;
    std::string mode_;
    double oracle_time_;
    double ve_time_;
    int no_facets_;

};

#endif /* defined(__mco__ap_benson_module__) */
