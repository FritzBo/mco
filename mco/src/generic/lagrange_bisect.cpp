//
//  lagrange_bisect.cpp
//  mco
//
//  Created by Fritz Bökler on 07.08.14.
//
//

#include <mco/generic/lagrange_bisect.h>

#include <iostream> // FIXME

using std::function;

namespace mco {

double LagrangeBisect::
find_multiplier(oracle_type oracle,
                unsigned dimension,
                const Point& bounds,
                Point& lambda) {
    
    for(unsigned i = 0; i < dimension; ++i) {
        lambda[i] = 1.0;
    }
    
    double value = 0, new_value = 0;
    bool done = false;
    for(unsigned i = 0; !done && i < 10; ++i) {
        for(unsigned k = 0; k < dimension; ++k) {
            new_value = iteration(k,
                              dimension,
                              lambda,
                              bounds,
                              oracle);
            
        }
        
        if(std::abs(value - new_value) < 1E-5) {
            done = true;
        }
        
        value = new_value;
    }
    
    return value + lambda * bounds;
}
    
double LagrangeBisect::iteration(unsigned objective,
                                unsigned dimension,
                                Point& lambda,
                                const Point& bounds,
                                oracle_type oracle) {

    // Probe l = infty
    Point temp_lambda = Point(0.0, dimension);
    temp_lambda[objective] = 1;
//    std::cout << temp_lambda << std::endl;
    Point point_right(dimension);
    double value_right = oracle(temp_lambda, point_right) - temp_lambda * bounds;
//    std::cout << point_right << std::endl;

    
    // Probe l = 0
    temp_lambda = lambda;
    temp_lambda[objective] = 0;
//    std::cout << temp_lambda << std::endl;
    Point point_left(dimension);
    double value_left = oracle(temp_lambda, point_left) - temp_lambda * temp_lambda;
//    std::cout << point_left << std::endl;
    
    while(value_left != value_right &&
          point_right[objective] != point_left[objective]) {
    
        // Compute new l
        double new_lambda = recompute_lambda(objective,
                                             lambda,
                                             point_left,
                                             point_right);
        
        if(std::abs(new_lambda - temp_lambda[objective]) < 1E-5) {
            break;
        }
        
        temp_lambda[objective] = new_lambda;
//        std::cout << temp_lambda << std::endl;
    
        // probe with new l
        Point new_point(dimension);
        double new_value = oracle(temp_lambda, new_point) - temp_lambda * bounds;
//        std::cout << new_point << std::endl;
        
        // replace new point and value
        if(new_point[objective] <= bounds[objective]) {
            
            point_right = new_point;
            value_right = new_value;
            value_left = temp_lambda * (point_left - bounds);
        } else {
            point_left = new_point;
            value_left = new_value;
            value_right = temp_lambda * (point_right - bounds);
        }
    }
    
    lambda = temp_lambda;
    return std::max(value_left, value_right);
}
    
double LagrangeBisect::recompute_lambda(unsigned objectve,
                                        const Point& lambda,
                                        const Point& point1,
                                        const Point& point2) {
    Point diff_point = point1 - point2;
    
    double x1 = 0;
    for(unsigned i = 0; i < point1.dimension(); ++i) {
        if(i != objectve) {
            x1 += lambda[i] * diff_point[i];
        }
    }
    
    return - x1/diff_point[objectve];
}
    
}