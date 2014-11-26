//
//  sp_lagrange_bisect.h
//  mco
//
//  Created by Fritz Bökler on 07.08.14.
//
//

#ifndef mco_sp_lagrange_bisect_h
#define mco_sp_lagrange_bisect_h

#include <functional>

#include <mco/basic/point.h>

namespace mco {

class LagrangeBisect {
    using oracle_type = std::function<double(const Point&, Point&)>;
public:
    double find_multiplier(oracle_type oracle,
                         unsigned dimension,
                         const Point& bounds,
                         Point& lambda);
    
private:
    double iteration(unsigned objective,
                   unsigned dimension,
                   Point& lambda,
                   const Point& bounds,
                   oracle_type oracle);
    
    double recompute_lambda(unsigned objective,
                            const Point& lambda,
                            const Point& point1,
                            const Point& point2);
    
};
    
}

#endif
