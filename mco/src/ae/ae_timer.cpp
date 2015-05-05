//
//  ae_timer.cpp
//  mco
//
//  Created by Fritz Bökler on 01.05.15.
//
//

#include <mco/ae/ae_timer.h>

using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

namespace mco {

std::ostream& operator<<(std::ostream& stream, const AeTimer& timer) {
    stream << "[";
    for(auto z : timer.durations()) {
        stream << z << ", ";
    }
    stream << "]";
    
    return stream;
}

}