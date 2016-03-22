//
//  ep_lc_approx_module.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__ep_bs_par_module__
#define __mco__ep_bs_par_module__

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EpBsParModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    EpBsParModule()
    :   AlgorithmModule(OGDF) { }

    virtual void perform(int argc, char** args);
    virtual ~EpBsParModule() {}
    
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
        
    void calculate_ideal_heuristic(const ogdf::Graph& graph,
                                   const ogdf::EdgeArray<mco::Point>& costs,
                                   unsigned dimension,
                                   const ogdf::node source,
                                   const ogdf::node target,
                                   bool directed,
                                   std::vector<ogdf::NodeArray<double>>& distances);

    static void first_phase(const ogdf::Graph* graph,
                            std::function<mco::Point*(ogdf::edge)> cost_function,
                            unsigned dimension,
                            const ogdf::node source,
                            const ogdf::node target,
                            std::function<void(ogdf::NodeArray<mco::Point*>&, ogdf::NodeArray<ogdf::edge>&)> callback);

};

#endif /* defined(__mco__ep_bs_module__) */
