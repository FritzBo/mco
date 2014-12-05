/*!
 \file
 \brief Contains code for comparing two Point objects by the lexicographic ordering.
 \date 2014-03-31
 \author Fritz Bökler

*/

#ifndef mco_lex_point_comparator_h
#define mco_lex_point_comparator_h

namespace mco {

// Class definitions

/*!
 Functor which compares to Point objects by the lexicographic order.

 \date 2014-03-31
 \author Fritz Bökler
 \ingroup basic
 */
class LexPointComparator {
public:
    /*!
     \param epsilon The epsilon value to use in comparisons.
     \param strict If true (default), the strict lexicographic ordering is used.
     */
    inline LexPointComparator(double epsilon = 0, bool strict = true)
    : epsilon_(epsilon), strict_(strict) {}
    
    inline bool operator()(const Point * p1,
                           const Point * p2) const noexcept;
    
    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;
    
    inline static bool is_lex_le(const Point& p1, const Point& p2, double epsilon) noexcept;
    
    inline static bool is_lex_leq(const Point& p1, const Point& p2, double epsilon) noexcept;
    
private:
    const double epsilon_;
    const bool strict_;
};

/*!
 \brief Functor which compares two Point objects by the lexicographic ordering where a component before the first is added which corresponds to the weighted sum of the objectives.
 
 \date 2014-03-31
 \author Fritz Bökler
 \ingroup basic
 */
class LexWeightedPointComparator {
public:
    /*!
     \param weighting A Point object such that the inner product of weighting and the given Points will form a new first component.
     \param epsilon The epsilon value to use in comparisons.
     \param strict If true (default), the strict lexicographic ordering is used.
     */
    inline LexWeightedPointComparator(Point weighting,
                                      double epsilon = 0,
                                      bool strict = true)
    :   epsilon_(epsilon),
        strict_(strict),
        weighting_(weighting) { }

    inline bool operator()(const Point * p1,
                           const Point * p2) const noexcept;

    inline bool operator()(const Point& p1,
                           const Point& p2) const noexcept;


private:
    const double epsilon_;
    const bool strict_;

    Point weighting_;
};

// Implementation of inline methods.

inline bool LexPointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {
    
    if(strict_) {
        return is_lex_le(p1, p2, epsilon_);
    } else {
        return is_lex_leq(p1, p2, epsilon_);
    }

}

inline bool LexPointComparator::
operator()(const Point * point1,
           const Point * point2) const noexcept {
    
    return operator()(*point1, *point2);
}

inline bool LexPointComparator::
is_lex_le(const Point &p1,
          const Point &p2,
          double epsilon) noexcept {
    
    assert(p1.dimension() == p2.dimension());
    
    unsigned dimension = p1.dimension();
    
    for(unsigned int i = 0; i < dimension; ++i) {
        if(p1[i] - p2[i] < - epsilon) {
            return true;
        } else if(p1[i] - p2[i] > epsilon) {
            return false;
        }
    }
    return false;
}

inline bool LexPointComparator::
is_lex_leq(const Point &p1,
           const Point &p2,
           double epsilon) noexcept {
    
    assert(p1.dimension() == p2.dimension());
    
    unsigned dimension = p1.dimension();
    
    for(unsigned int i = 0; i < dimension; ++i) {
        if(p1[i] - p2[i] < - epsilon) {
            return true;
        } else if(p1[i] - p2[i] > epsilon) {
            return false;
        }
    }
    return true;
}

inline bool LexWeightedPointComparator::
operator()(const Point& p1,
           const Point& p2) const noexcept {
    
    double weight1 = p1 * weighting_;
    double weight2 = p2 * weighting_;
    
    if(weight1 < weight2) {
        return true;
    } else if(weight2 < weight1) {
        return false;
    }
    
    if(strict_) {
        return LexPointComparator::is_lex_le(p1, p2, epsilon_);
    } else {
        return LexPointComparator::is_lex_leq(p1, p2, epsilon_);
    }
    
}

inline bool LexWeightedPointComparator::
operator()(const Point* p1,
           const Point* p2) const noexcept {
    
    return operator()(*p1, *p2);
}

}

#endif