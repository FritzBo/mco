//
//  cdd_files.h
//  mco
//
//  Created by Fritz Bökler on 26.11.14.
//
//

#ifndef mco_cdd_files_h
#define mco_cdd_files_h

#include <ostream>

#include <mco/basic/point.h>

namespace mco {

class CddFiles {
public:
    template<typename InputIterator>
    static void write_v_rep(std::ostream& stream,
                            unsigned size,
                            unsigned dimension,
                            InputIterator begin,
                            InputIterator end);

    template<typename InputIterator>
    static void write_h_rep(std::ostream& stream,
                            unsigned size,
                            unsigned dimension,
                            InputIterator begin,
                            InputIterator end);
};

template<typename InputIterator>
void CddFiles::write_v_rep(std::ostream& stream,
                           unsigned size,
                           unsigned dimension,
                           InputIterator begin,
                           InputIterator end) {

    stream << "V-representation" << endl;
    stream << "begin" << endl;
    stream << (size + dimension) << "\t" << (dimension + 1) << "\t" << "real" << endl;

    // Printing extreme rays (Pareto-Cone)
    for(unsigned i = 0; i < dimension; ++i) {
        stream << "0";
        for(unsigned j = 0; j < dimension; ++j) {
            stream << "\t";
            if(i != j) {
                stream << "0";
            } else {
                stream << "1";
            }
        }
        stream << endl;
    }

    // Printing extreme points
    auto it = begin;
    while(it != end) {

        const Point& p = it->second;

        stream << "1";

        for(unsigned i = 0; i < dimension; ++i) {
            stream << "\t" << p[i];
        }
        stream << endl;

        ++it;
    }

    stream << "end" << endl;

}

template<typename InputIterator>
void CddFiles::write_h_rep(std::ostream& stream,
                           unsigned size,
                           unsigned dimension,
                           InputIterator begin,
                           InputIterator end) {

    stream << "H-representation" << endl;
    stream << "begin" << endl;
    stream << size << "\t" << dimension << "real" << endl;

    auto it = begin;
    while(it != end) {

        const Point& p = it->second;

        stream << (- p[dimension - 1]);

        for(unsigned i = 0; i < dimension - 1; ++i) {
            stream << "\t" << p[i];
        }

        stream << endl;

        ++it;
    }

    stream << "end" << endl;
}

}

#endif
