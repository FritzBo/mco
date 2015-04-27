//
//  ep_tsaggouris_module.h
//  mco
//
//  Created by Fritz Bökler on 14.04.15.
//
//

#ifndef __mco__ep_tsaggouris_module__
#define __mco__ep_tsaggouris_module__

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EpTsaggourisModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~EpTsaggourisModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;
    
    void parse_ideal_bounds(const TCLAP::MultiArg<std::string>& argument,
                            unsigned dimension,
                            std::function<double(ogdf::node, unsigned)> heuristic,
                            const ogdf::node source,
                            mco::Point& bounds);

    void parse_absolute_bounds(const TCLAP::MultiArg<std::string>& argument,
                            unsigned dimension,
                            mco::Point& bounds);
    
    void parse_epsilon(const TCLAP::MultiArg<string>& epsilon_argument,
                       unsigned dimension,
                       mco::Point& epsilon);
    
    void calculate_ideal_heuristic(const ogdf::Graph& graph,
                                   const ogdf::EdgeArray<mco::Point>& costs,
                                   unsigned dimension,
                                   const ogdf::node source,
                                   const ogdf::node target,
                                   bool directed,
                                   std::vector<ogdf::NodeArray<double>>& distances);
    
    
};

#endif /* defined(__mco__ep_lc_approx_module__) */