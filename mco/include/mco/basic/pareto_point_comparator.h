/*!
 \file
 \brief Contains a functor to compare two Point objects using the Pareto-dominance relation.
 \author Fritz Bökler
 \date 2014-03-31
 */

#ifndef mco_pareto_point_comparator_h
#define mco_pareto_point_comparator_h

namespace mco {

/*!
 \brief A functor to compare two Point objects using the Pareto-dominace relation.
 \date 2014-03-31
 \author Fritz Bökler
 \ingroup basic
 */
class ParetoDominationPointComparator {
public:
    /*!
     \param epsilon The epsilon value used in the comparision. Defaults to 0.
     */
    inline ParetoDominationPointComparator(double epsilon = 0)
    : epsilon_(epsilon) { }

    //! Returns \c true iff \p p1 dominates \p p2.
    inline bool operator()(const Point * p1,
                           const Point * p2) const noexcept;

    //! Returns \c true iff \p p1 dominates \p p2.
    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;

    /*! \brief Static method which is used by operator(). Returns \c true iff \p p1 dominates \p p2 using the value \p epsilon for inexact comparisons.
     */
    inline static bool dominates(const Point& p1,
                                 const Point& p2,
                                 double epsilon) noexcept;
    
private:
    const double epsilon_;
};

inline bool ParetoDominationPointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {
    
    return dominates(p1, p2, epsilon_);
}

inline bool ParetoDominationPointComparator::
operator()(const Point * p1,
           const Point * p2) const noexcept {
    
    return dominates(*p1, *p2, epsilon_);
}

inline bool ParetoDominationPointComparator::
dominates(const Point &p1,
          const Point &p2,
          double epsilon) noexcept {
    
    assert(p1.dimension() == p2.dimension());
    
    unsigned dimension = p1.dimension();
    
    bool equal = true;
    
    for(unsigned int i = 0; i < dimension; i++) {
        
        if(p1[i] - p2[i] > epsilon) {
            return false;
        }
        
        if(abs(p1[i] - p2[i]) > epsilon) {
            equal = false;
        }
    }
    
    return !equal;
}

}

#endif
