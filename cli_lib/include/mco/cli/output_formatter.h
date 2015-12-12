//
//  output_formatter.h
//  mco
//
//  Created by Fritz Bökler on 04.08.14.
//
//

#ifndef __mco__output_formatter__
#define __mco__output_formatter__

#include <iterator>
#include <iostream>

#include <tclap/CmdLine.h>

namespace mco {
    
template<typename S>
class CliOutputFormatter {
    
    using solution_type = S;

public:
    CliOutputFormatter(TCLAP::CmdLine* cmd_line)
    :   cmd_line_(cmd_line) { }
    
    inline bool initialize_cmd_line();
    
    template<typename InputIterator>
    inline bool print_output(std::ostream& ostream,
                             InputIterator begin,
                             InputIterator end);
    
    inline ~CliOutputFormatter();
    
private:
    TCLAP::CmdLine* cmd_line_;
    
    TCLAP::SwitchArg* print_frontier_arg_;
    TCLAP::SwitchArg* print_solutions_arg_;
    TCLAP::SwitchArg* force_print_all_arg_;
    TCLAP::SwitchArg* print_count_arg_;
    TCLAP::SwitchArg* print_timing_arg_;
    TCLAP::SwitchArg* print_statistics_;
};

template<typename S>
CliOutputFormatter<S>::~CliOutputFormatter() {
    delete print_frontier_arg_;
    delete print_solutions_arg_;
    delete force_print_all_arg_;
    delete print_count_arg_;
    delete print_timing_arg_;
    delete print_statistics_;
}
    
template<typename S>
inline bool CliOutputFormatter<S>::
initialize_cmd_line() {
    print_frontier_arg_ =
        new TCLAP::SwitchArg("F",
                             "print-frontier",
                             "Prints the found frontier if "
                             "the problem was feasible.",
                             false);
    
    print_solutions_arg_ =
        new TCLAP::SwitchArg("S",
                             "print-solutions",
                             "Print the found solutions if "
                             "the problem was feasible.",
                             false);
    
    force_print_all_arg_ =
        new TCLAP::SwitchArg("A",
                             "force-all",
                             "Prints all points/solutions, no matter how many",
                             false);
    
    print_count_arg_ =
        new TCLAP::SwitchArg("C",
                             "count",
                             "Prints the size of the Pareto-front",
                             false);
    
    print_timing_arg_ =
        new TCLAP::SwitchArg("T",
                             "timing",
                             "Prints timing information",
                             false);
    
    print_statistics_ =
        new TCLAP::SwitchArg("I",
                             "statistics",
                             "Prints statistics of the algorithm",
                             false);
    
    cmd_line_->add(print_frontier_arg_);
    cmd_line_->add(print_solutions_arg_);
    cmd_line_->add(force_print_all_arg_);
    cmd_line_->add(print_count_arg_);
    cmd_line_->add(print_timing_arg_);
    cmd_line_->add(print_statistics_);
    
    return true;
}
    
template<typename S>
template<typename ConstInputIterator>
inline bool CliOutputFormatter<S>::
print_output(std::ostream& ostream,
             ConstInputIterator begin,
             ConstInputIterator end) {
    
    bool print_frontier     = print_frontier_arg_->getValue();
    bool print_solutions    = print_solutions_arg_->getValue();
    bool print_count        = print_count_arg_->getValue();
    bool force_print_all    = force_print_all_arg_->getValue();
    bool print_timing       = print_timing_arg_->getValue();
    bool print_statistics   = print_statistics_->getValue();
    
    if(print_frontier || print_solutions || print_count) {
        
        unsigned frontier_size = std::distance(begin, end);
        
        ostream << frontier_size << std::endl;
        
        if(print_frontier || print_solutions) {
            
            int count = 0;
            auto solution_it = begin;
            while(solution_it != end && (count < 25 || force_print_all)) {
                auto solution = *solution_it;
                
                if(print_frontier) {
                    auto point_it = solution.second.cbegin();
                    while(point_it != solution.second.cend()) {
                        ostream << *point_it++ << ", ";
                    }
                }
                if(print_solutions) {
                    for(auto edge : solution.first) {
                        ostream << edge->source() << ", ";
                    }
                    ostream << (*solution.first.rbegin())->target() << ", ";
                }
                std::cout << std::endl;
                solution_it++;
                count++;
            }
            
        }
    }
    
    if(print_statistics) {
        std::cout << "Some statistics information" << std::endl;
    }
    
    if(print_timing) {
        std::cout << "Some timining information" << std::endl;
    }

    
    return true;
    
}
    
}

#endif /* defined(__mco__output_formatter__) */
