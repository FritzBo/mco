//
//  parse_util.h
//  mco
//
//  Created by Fritz Bökler on 06.08.14.
//
//

#ifndef mco_parse_util_h
#define mco_parse_util_h

#include <functional>

#include <tclap/CmdLine.h>

#include <ogdf/basic/Graph.h>

namespace mco {
    
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

inline
bool parse_ideal_bounds(const TCLAP::MultiArg<string>& argument,
                        unsigned dimension,
                        std::function<double(ogdf::node, unsigned)> heuristic,
                        const ogdf::node source,
                        Point& bound) {
    
    auto bounds_it = argument.begin();
    while(bounds_it != argument.end()) {
        vector<string> tokens;
        tokenize(*bounds_it, tokens, ":");
        
        if(tokens.size() != 2) {
            return false;
        }
        
        unsigned objective_function = stoul(tokens[0]);
        double factor = stod(tokens[1]);
        
        if(objective_function > dimension) {
            return false;
        }
        
        if(objective_function == 0) {
            return false;
        }
        
        bound[objective_function - 1] = factor * heuristic(source, objective_function - 1);
        
        bounds_it++;
    }
    
    return true;
}
    
}

#endif
