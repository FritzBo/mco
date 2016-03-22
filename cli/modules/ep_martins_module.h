//
//  benson_module.h
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#ifndef __mco__martins_module__
#define __mco__martins_module__

#include <list>
#include <string>

#include <ogdf/basic/Graph.h>

#include "../basic/modules.h"

class EpMartinsModule : public AlgorithmModule<std::list<ogdf::edge>> {
    
public:
    EpMartinsModule()
    :   AlgorithmModule(OGDF) { }

    virtual void perform(int argc, char** args);
    virtual ~EpMartinsModule() {}
    
    virtual const std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>>& solutions();
    virtual std::string statistics();
    
private:
    
    void calculate_ideal_heuristic(const ogdf::Graph& graph,
                                   const ogdf::EdgeArray<mco::Point>& costs,
                                   unsigned dimension,
                                   const ogdf::node source,
                                   const ogdf::node target,
                                   bool directed,
                                   std::vector<ogdf::NodeArray<double>>& distances);
    
    void parse_ideal_bounds(const TCLAP::MultiArg<std::string>& argument,
                            unsigned dimension,
                            std::function<double(ogdf::node, unsigned)> heuristic,
                            const ogdf::node source,
                            std::list<mco::Point> & bounds);
    
    void parse_fractional_bounds(const TCLAP::MultiArg<std::string>& argument,
                                 unsigned dimension,
                                 mco::Point& bounds);
    
    void parse_bundling_bound(TCLAP::ValueArg<std::string>& bundling_arg,
                              unsigned dimension,
                              std::function<double(ogdf::node, unsigned)> ideal_heuristic,
                              const mco::Point& ideal_bound,
                              mco::Point& bundling_bound);

    
    static void first_phase(const ogdf::Graph* graph,
                            std::function<mco::Point*(ogdf::edge)> cost_function,
                            unsigned dimension,
                            const ogdf::node source,
                            const ogdf::node target,
                            double epsilon,
                            std::function<void(ogdf::NodeArray<mco::Point*>&, ogdf::NodeArray<ogdf::edge>&)> callback);
    
    std::list<std::pair<const std::list<ogdf::edge>, const mco::Point>> solutions_;

    
};

template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", bool trimEmpty = false)
{

    using ValueType = typename ContainerT::value_type;
    using SizeType = typename ValueType::size_type;
    
    std::string::size_type pos, lastPos = 0;
    while(true)
    {
        pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos)
        {
            pos = str.length();
            
            if(pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos,
                                           (SizeType) pos - lastPos ));
            
            break;
        }
        else
        {
            if(pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos,
                                           (SizeType)pos - lastPos ));
        }
        
        lastPos = pos + 1;
    }
};

#endif /* defined(__mco__martins_module__) */
