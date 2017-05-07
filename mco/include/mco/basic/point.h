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
#include <vector>

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
#ifdef NDEBUG
    virtual ~Point() noexcept { delete[] values_; }
#endif

    //! Creates a well defined Point object with no components.
#ifdef NDEBUG
    inline Point() : dimension_(0), values_(nullptr) { }
#else
    inline Point() : dimension_(0), values_(0) { }
#endif

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
#ifdef NDEBUG
    const double* cbegin() const { return values_; }
#else
    std::vector<double>::const_iterator cbegin() const { return values_.cbegin(); }
#endif

    //! Returns a const_iterator to end.
#ifdef NDEBUG
    const double* cend() const { return values_ + dimension_; }
#else
    std::vector<double>::const_iterator cend() const { return values_.cend(); }
#endif

    //! Returns an interator to beginning.
#ifdef NDEBUG
    double* begin() { return values_; }
#else
    std::vector<double>::iterator begin() { return values_.begin(); }
#endif

    //! Returns an iterator to end.
#ifdef NDEBUG
    double* end() { return values_ + dimension_; }
#else
    std::vector<double>::iterator end() { return values_.end(); }
#endif

    friend bool operator==(Point const &, Point const &);
    friend bool operator!=(Point const &, Point const &);
	friend std::ostream & operator<<(std::ostream &, const Point &);
    friend void swap(Point& p1, Point& p2);

    
protected:
    unsigned int dimension_ = 0; //! The number of components of this Point object.

#ifdef NDEBUG
    double * values_ = nullptr; //! A pointer pointing to the values of this point object.
#else
    std::vector<double> values_;
#endif
};

inline Point::Point(unsigned dimension)
: dimension_(dimension)
#ifndef NDEBUG
    , values_(dimension)
#endif
{
#ifdef NDEBUG
    values_ = dimension_ ? new double[dimension] : nullptr;
#endif
    for(unsigned int i = 0; i < dimension; ++i)
        values_[i] = 0;
}
    
inline Point::Point(const double *values, unsigned int dimension)
: dimension_(dimension)
#ifndef NDEBUG
    , values_(dimension)
#endif
{
#ifdef NDEBUG
    values_ = dimension_ ? new double[dimension] : nullptr;
#endif
    std::copy(values, values + dimension, begin());
}
    
inline Point::Point(double value, unsigned int dimension)
: dimension_(dimension)
#ifndef NDEBUG
    , values_(dimension)
#endif
{
#ifdef NDEBUG
    values_ = dimension_ ? new double[dimension] : nullptr;
#endif
    for(unsigned i = 0; i < dimension; ++i) {
        values_[i] = value;
    }
}
    
inline Point::Point(std::initializer_list<double> values)
: dimension_(values.size())
#ifndef NDEBUG
    , values_(values.size())
#endif
{
#ifdef NDEBUG
    values_ = dimension_ ? new double[values.size()] : nullptr;
#endif
    std::copy(values.begin(), values.end(), begin());
}

inline Point::Point(const Point &p) : dimension_(p.dimension_) {
#ifdef NDEBUG
    values_ = dimension_ ? new double[dimension_] : nullptr;
#else
    values_.resize(dimension());
#endif
    std::copy(p.cbegin(), p.cend(), begin());
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
    
    return that;
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

//! Two point objects are equal, if they point to the same array of points.
inline bool operator==(Point const & p1, Point const & p2)
{
    return p1.values_ == p2.values_;
}

//! Two point objects are equal, if they point to the same array of points.
inline bool operator!=(Point const & p1, Point const & p2)
{
    return !(p1 == p2);
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

class PointSet {
public:
    PointSet() = delete;

    //! Creates a new point set of points with \p dimension components.
    explicit inline PointSet(unsigned dimension)
    :   dimension_(dimension),
        point_array_() { }

    //! Creates a new point set of size \p size of points with \p dimension components
    inline PointSet(unsigned dimension, unsigned size)
    :   dimension_(dimension),
        point_array_(size * dimension, 0.0) { }

    //! Removes a point from the set
    inline void remove(unsigned index);

    //! Create a new point with zero components.
    inline void new_point();

private:
    unsigned dimension_;
    std::vector<double> point_array_;
};

void PointSet::remove(unsigned index)
{
    unsigned copyend = point_array_.size() - 1;
    unsigned deletend = index + dimension_;

    for(unsigned i = 0; i < dimension_; ++i)
    {
        point_array_[deletend] = point_array_[copyend];
        copyend -= 1;
        deletend -= 1;
    }

    point_array_.resize(point_array_.size() - dimension_);
}

void PointSet::new_point()
{
    point_array_.resize(point_array_.size() + dimension_, 0.0);
}

class PointReference : public Point {

    //! Destructor removes point from PointSet
    virtual ~PointReference() noexcept
    {
        reference_set_.remove(index_);
    }

    //! Creates a well defined Point object with zeros.
    inline PointReference(PointSet& ref_set);

    //! Creates a Point object with \p dimension components consisting of the first \p dimension values of \p values.
    inline PointReference(const double *values, unsigned int dimension);

    //! Creates a Point object with \p dimension components which are all set to \p value.
    inline PointReference(double value, unsigned int dimension);

    //! Creates a Point object using all values from the initialzer list values.
    inline PointReference(std::initializer_list<double> values);

    //! Copy constructor
    inline PointReference(const Point &p);

    //! Move constructor (employing swap)
    inline PointReference(Point&& that) noexcept;

    //! Copy assignment
    inline PointReference & operator=(const Point& that) noexcept;

    //! Move assignment (employing swap)
    inline PointReference & operator=(Point&& that) noexcept;
    
protected:
    PointSet& reference_set_;
    unsigned index_;
};

PointReference::PointReference(PointSet& ref_set)
:   Point(),
    reference_set_(ref_set)
{
    reference_set_.new_point();
}


} // namespace mco

#include "lex_point_comparator.h"
#include "equality_point_comparator.h"
#include "componentwise_point_comparator.h"
#include "pareto_point_comparator.h"


#endif /* MCO_POINT_H_ */
