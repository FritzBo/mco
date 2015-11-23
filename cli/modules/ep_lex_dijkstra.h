//
//  benson_module.h
//  mco
//
//  Created by Fritz Bökler on 01.09.14.
//
//

#ifndef __mco__ld_module__
#define __mco__ld_module__

#include <list>

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EpLDModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~EpLDModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;

    mco::Point parse_weight(TCLAP::ValueArg<string>& argument,
                            unsigned dimension);

};

#endif /* defined(__mco__ld_module__) */
