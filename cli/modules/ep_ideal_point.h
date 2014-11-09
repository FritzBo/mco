//
//  benson_module.h
//  mco
//
//  Created by Fritz Bökler on 01.09.14.
//
//

#ifndef __mco__ideal_module__
#define __mco__ideal_module__

#include <list>

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EpIdealModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~EpIdealModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;

};

#endif /* defined(__mco__martins_module__) */