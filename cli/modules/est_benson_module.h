//
//  benson_module.h
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#ifndef __mco__est_benson_module__
#define __mco__est_benson_module__

#include <list>

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EstBensonModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~EstBensonModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;

    unsigned dimension;
    unsigned nodes_;
    unsigned edges_;


};

#endif /* defined(__mco__est_benson_module__) */
