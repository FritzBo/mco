/*!
 \file
 \brief Contains a functor for the canonical partial order on R^d.
 \date 2014-03-31
 \author Fritz Bökler
 */

#ifndef mco_componentwise_point_comparator_h
#define mco_componentwise_point_comparator_h

namespace mco {

// Class Definition

/*!
 \brief A functor to compare two point objects using the canonical partial order on R^d.
 \author Fritz Bökler
 \date 2014-03-31
 \ingroup basic
 */
class ComponentwisePointComparator {
public:
    /*!
     \param epsilon The epsilon value to be used in comparison. Defaults to 0.
     \param strict If true (default), the strict partial order is used.
     */
    inline ComponentwisePointComparator(double epsilon = 0, bool strict = true)
    : epsilon_(epsilon), strict_(strict) { }
    
    inline bool operator()(const Point * point1,
                           const Point * point2) const noexcept;
    
    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;
    
    inline static bool is_le(const Point& p1,
                             const Point& p2,
                             double epsilon) noexcept;
    
    inline static bool is_leq(const Point& p1,
                              const Point& p2,
                              double epsilon) noexcept;
    
private:
    const double epsilon_;
    const bool strict_;
};

// Implementation of inline functions

inline bool ComponentwisePointComparator::
operator()(const Point* p1,
           const Point* p2) const noexcept {
    
    if(strict_) {
        return is_le(*p1, *p2, epsilon_);
    } else {
        return is_leq(*p1, *p2, epsilon_);
    }
}

inline bool ComponentwisePointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {
    
    if(strict_) {
        return is_le(p1, p2, epsilon_);
    } else {
        return is_leq(p1, p2, epsilon_);
    }
}

inline bool ComponentwisePointComparator::
is_le(const Point &p1,
      const Point &p2,
      double epsilon) noexcept {
    
    assert(p1.dimension() == p2.dimension());
    
    unsigned dimension = p1.dimension();
    
    for(unsigned int i = 0; i < dimension; ++i) {
        if(p1[i] - p2[i] >= - epsilon) {
            return false;
        }
    }
    
    return true;
}


inline bool ComponentwisePointComparator::
is_leq(const Point &p1,
       const Point &p2,
       double epsilon) noexcept {
    
    assert(p1.dimension() == p2.dimension());
    
    unsigned dimension = p1.dimension();
    
    for(unsigned int i = 0; i < dimension; ++i) {
        if(p1[i] - p2[i] > epsilon) {
            return false;
        }
    }
    
    return true;
}

}

#endif
