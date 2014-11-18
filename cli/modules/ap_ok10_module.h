//
//  benson_module.h
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#ifndef __mco__ap_ok10_module__
#define __mco__ap_ok10_module__

#include <list>

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class ApOk10Module : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~ApOk10Module() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;

    
};

#endif /* defined(__mco__ap_ok10_module__) */
