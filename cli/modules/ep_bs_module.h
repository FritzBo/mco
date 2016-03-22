//
//  ep_lc_approx_module.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__ep_bs_module__
#define __mco__ep_bs_module__

#include <mco/basic/forward_star.h>

#include "../basic/modules.h"

class EpBsModule : public AlgorithmModule<std::list<mco::node>> {
    
public:
    EpBsModule()
    :   AlgorithmModule(FORWARDSTAR) { }

    virtual void perform(int argc, char** args);
    virtual ~EpBsModule() {}
    
    virtual const std::list<std::pair<const std::list<mco::node>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:

    // stats:
    unsigned        num_nodes_                          = 0;
    unsigned        num_edges_                          = 0;
    unsigned        num_objectives_                     = 0;
    unsigned long   label_compares_                     = 0;
    unsigned        deleted_tree_labels_                = 0;
    unsigned        recursive_deletions_                = 0;
    unsigned        arc_pushes_                         = 0;
    unsigned        touched_recursively_deleted_label_  = 0;
    unsigned long   deleted_labels_                     = 0;

    double          solution_time_                      = 0;
    
    std::list<std::pair<const std::list<mco::node>, const mco::Point>> solutions_;

    bool label_select_;
    
//    void parse_ideal_bounds(const TCLAP::MultiArg<std::string>& argument,
//                            unsigned dimension,
//                            std::function<double(ogdf::node, unsigned)> heuristic,
//                            const ogdf::node source,
//                            mco::Point& bounds);

    void parse_absolute_bounds(const TCLAP::MultiArg<std::string>& argument,
                            unsigned dimension,
                            mco::Point& bounds);
        
//    void calculate_ideal_heuristic(const ogdf::Graph& graph,
//                                   const ogdf::EdgeArray<mco::Point>& costs,
//                                   unsigned dimension,
//                                   const ogdf::node source,
//                                   const ogdf::node target,
//                                   bool directed,
//                                   std::vector<ogdf::NodeArray<double>>& distances);

//    static void first_phase(const ogdf::Graph* graph,
//                            std::function<mco::Point*(ogdf::edge)> cost_function,
//                            unsigned dimension,
//                            const ogdf::node source,
//                            const ogdf::node target,
//                            std::function<void(ogdf::NodeArray<mco::Point*>&, ogdf::NodeArray<ogdf::edge>&)> callback);

};

#endif /* defined(__mco__ep_bs_module__) */
