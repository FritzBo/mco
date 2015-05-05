//
//  ae_timer.h
//  mco
//
//  Created by Fritz Bökler on 01.05.15.
//
//

#ifndef __mco__ae_timer__
#define __mco__ae_timer__

#include <chrono>
#include <vector>
#include <ostream>

namespace mco {

class AeTimer {
public:

    inline void start() {
        anchor_ = std::chrono::steady_clock::now();
    }

    inline void stop() {
        log_time();
        anchor_ = std::chrono::steady_clock::now();
    }

    inline void log_time() {
        using std::chrono::steady_clock;
        using std::chrono::duration_cast;
        using std::chrono::duration;

        auto now = steady_clock::now();
        durations_.push_back(duration_cast<duration<double>>(now - anchor_).count());
    }

    inline const std::vector<double> & durations() const {
        return durations_;
    }

    inline double last() const {
        return durations_.back();
    }

    inline double sum() const {
        double sum = 0;
        for(auto z : durations_) {
            sum += z;
        }
        return sum;
    }
    
private:
    std::vector<double> durations_;
    std::chrono::steady_clock::time_point anchor_;
};

std::ostream& operator<<(std::ostream& stream, const AeTimer& timer);

}

#endif /* defined(__mco__ae_timer__) */
