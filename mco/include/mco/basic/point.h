#pragma once
/*! \file
 \brief Holds all information for the mco::Point class and additional functions usuful in connection with the mco::Point class.
 \date 2013-03-15
 \author Fritz Bökler
 */

#ifndef MCO_POINT_H_
#define MCO_POINT_H_

#include <ostream>
#include <stdexcept>
#include <cassert>
#include <cmath>

namespace mco {

//! Class to model points in R^d.
/*!
 Exceptions are omitted for performance reasons. All methods are inlined to also allow for better performance.

 \author Fritz Bökler
 \ingroup basic
 \date 2013-03-15
 \todo Memory Pool allocator for the chunks of IEEE data.
 */
class Point {
public:

    //! Destructor deletes all allocated values.
    virtual ~Point() noexcept { delete[] values_; }

    //! Creates a well defined Point object with no components.
    inline Point() : dimension_(0), values_(nullptr) { }

    //! Creates a well defined Point object with \p dimension components. The components will be initialized with 0.0.
    inline explicit Point(unsigned int dimension);

    //! Creates a Point object with \p dimension components consisting of the first \p dimension values of \p values.
	inline Point(const double *values, unsigned int dimension);

    //! Creates a Point object with \p dimension components which are all set to \p value.
    inline Point(double value, unsigned int dimension);

    //! Creates a Point object using all values from the initialzer list values.
    inline Point(std::initializer_list<double> values);

    //! Copy constructor
	inline Point(const Point &p);
    
    //! Move constructor (employing swap)
    inline Point(Point&& that) noexcept;
    
    //! Copy assignment
	inline Point & operator=(const Point& that) noexcept;
    
    //! Move assignment (employing swap)
    inline Point & operator=(Point&& that) noexcept;

    //! Returns the number of components the Point object has.
    unsigned dimension() const noexcept { return dimension_; }

    //! Adds the Point object \p that to this object.
    /*!
     In debug mode, an assertion is fired if both objects have different number of components.

     \warning In release mode the compatibility of the numbers of components will not be checked!
     */
    inline Point & operator+=(const Point &that) noexcept;

    //! Substracts \p that from the point object.
    /*!
     In debug mode, an assertion is fired if both objects have different number of components.

     \warning In release mode the compatibility of the numbers of components will not be checked!     */
    inline Point & operator-=(const Point &that) noexcept;

    //! Conducts a multiplication of this Point object with the scalar \p d.
    inline Point & operator*=(double d) noexcept;

    //! Returns a new Point object which represents the negative of this Point object.
    inline Point operator-() const &;

    //! Sames as (*this) *= -1;
    inline Point operator-() &&;

    //! Adds this Point object and \p that.
    /*!
     This version of addition works if that is an rvalue or called by value.

     In debug mode, an assertion is fired if both objects have different number of components.

     \warning In release mode the compatibility of the numbers of components will not be checked!     */
    inline Point operator+(Point that) const & noexcept;

    //! Adds this Point object and \p that.
    /*!
     This version of addition works if this is an rvalue and that is an lvalue.

     In debug mode, an assertion is fired if both objects have different number of components.

     \warning In release mode the compatibility of the numbers of components will not be checked!
     */
    inline Point operator+(const Point& that) && noexcept;

    //! Subtracts \p that from this.
    /*!
     This version of subtract works if that is an rvalue or called by value.

     In debug mode, an assertion is fired if both objects have different number of components.

     \warning In release mode the compatibility of the numbers of components will not be checked! 
     */
    inline Point operator-(Point that) const & noexcept;

    //! Subtracts \p that from this.
    /*!
     This version of subtract works if this is an rvalue and that is an lvalue.

     In debug mode, an assertion is fired if both objects have different number of components.

     \warning In release mode the compatibility of the numbers of components will not be checked!
     */
    inline Point operator-(const Point& that) && noexcept;
    
    //! Inner product of this and \p that
    inline double operator*(const Point &point) const noexcept;

    //! Returns a reference to component \p index.
	double& operator[](unsigned int index) { return values_[index]; }

    //! Returns a const reference to component \p index.
	const double& operator[](unsigned int index) const { return values_[index]; }

    //! Returns a const_iterator to beginning.
    const double* cbegin() const { return values_; }

    //! Returns a const_iterator to end.
    const double* cend() const { return values_ + dimension_; }

    //! Returns an interator to beginning.
    double* begin() { return values_; }

    //! Returns an iterator to end.
    double* end() { return values_ + dimension_; }

	friend std::ostream & operator<<(std::ostream &, const Point &);
    friend void swap(Point& p1, Point& p2);

    
private:
    unsigned int dimension_ = 0; //! The number of components of this Point object.

    double * values_ = nullptr; //! A pointer pointing to the values of this point object.
};

inline Point::Point(unsigned dimension)
: dimension_(dimension) {
    values_ = dimension_ ? new double[dimension] : nullptr;
    for(unsigned int i = 0; i < dimension; ++i)
        values_[i] = 0;
}
    
inline Point::Point(const double *values, unsigned int dimension)
: dimension_(dimension) {
    values_ = dimension_ ? new double[dimension] : nullptr;
    std::copy(values, values + dimension, values_);
}
    
inline Point::Point(double value, unsigned int dimension)
: dimension_(dimension) {
    values_ = dimension_ ? new double[dimension] : nullptr;
    for(unsigned i = 0; i < dimension; ++i) {
        values_[i] = value;
    }
}
    
inline Point::Point(std::initializer_list<double> values)
: dimension_(values.size()) {
    values_ = dimension_ ? new double[values.size()] : nullptr;
    std::copy(values.begin(), values.end(), values_);
}

inline Point::Point(const Point &p) : dimension_(p.dimension_) {
    values_ = dimension_ ? new double[dimension_] : nullptr;
    std::copy(p.values_, p.values_ + dimension_, values_);
}
    
inline Point::Point(Point&& that) noexcept : Point() {
    swap(*this, that);
}
    
inline Point& Point::operator=(const Point& that) noexcept {
    
    if(that.dimension_ <= dimension_) {
        dimension_ = that.dimension_;
        std::copy(that.cbegin(), that.cend(), this->begin());
    } else {
        Point new_point(that);
        swap(*this, new_point);
    }
    
    return *this;
}
    
inline Point& Point::operator=(Point&& that) noexcept {
    swap(*this, that);
    return *this;
}

    
inline Point & Point::operator+=(const Point &p) noexcept {
    assert(dimension_ == p.dimension_);
    
    for(unsigned int i = 0; i < dimension_; ++i)
        values_[i] += p[i];
        
    return *this;
}

inline Point & Point::operator-=(const Point &p) noexcept {
    assert(dimension_ == p.dimension_);
    
    for(unsigned int i = 0; i < dimension_; ++i)
        values_[i] -= p[i];
        
    return *this;
}
    
inline Point & Point::operator*=(double d) noexcept {
    for(unsigned int i = 0; i < dimension_; ++i)
        values_[i] *= d;
        
    return *this;
}
    
inline Point Point::operator-() const & {
    Point result(*this);
    result *= -1;
    return result;
}
    
inline Point Point::operator-() && {
    *this *= -1;
    return std::move(*this);
}
    
inline Point Point::operator+(Point that) const & noexcept {
    assert(dimension_ == that.dimension_);
    return std::move(that += *this);
}

inline Point Point::operator+(const Point& that) && noexcept {
    assert(dimension_ == that.dimension_);
    return std::move(*this += that);
}
    
inline Point Point::operator-(Point that) const & noexcept {
    assert(dimension_ == that.dimension_);
    
    for(unsigned i = 0; i < dimension_; ++i) {
        that[i] = this->operator[](i) - that[i];
    }
    
    return std::move(that);
}
    
inline Point Point::operator-(const Point& that) && noexcept {
    assert(dimension_ == that.dimension_);
    return std::move(*this -= that);
}
    
inline double Point::operator*(const Point &point) const noexcept {
    assert(dimension_ == point.dimension_);
    
    double sum = 0;
    unsigned i;
    for(i = 1; i < dimension_; i += 2) {
        sum += (*this)[i - 1] * point[i - 1];
        sum += (*this)[i] * point[i];
    }
    
    return i == dimension_ + 1 ? sum : sum + (*this)[dimension_ - 1] * point[dimension_-1];
}

//! Multiplies the Point object \p p and the scalar value \p d.
inline Point operator*(Point p, double d) {
    return std::move(p *= d);
}

//! Multiplies the Point object \p p and the scalar value \p d.
inline Point operator*(double d, Point p) {
    return std::move(p *= d);
}

//! Writes a Point object to an output stream.
inline std::ostream & operator<<(std::ostream &os, const Point &point) {
    os << "(";
    
    if(point.dimension_ != 0) {
        for(unsigned int i = 0; i < point.dimension_ - 1; ++i) {
            os << point[i] << ", ";
        }
        
        os << point[point.dimension_ - 1];
    }
    
    os << ")";
    
    return os;
}

//! Swaps the value pointers of two Point objectives without copying.
inline void swap(Point& p1, Point& p2) {
    using std::swap;
    
    swap(p1.values_,    p2.values_);
    swap(p1.dimension_, p2.dimension_);
}

#include "lex_point_comparator"
#include "equality_point_comparator"
#include "componentwise_point_comparator"
#include "pareto_point_comparator"
    
} // namespace mco


#endif /* MCO_POINT_H_ */
