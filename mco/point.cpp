/*
 * point.cpp
 *
 *  Created on: 15.03.2013
 *      Author: fritz
 */

#include <mco/point.h>

#include <iostream>
#include <cstring>
#include <stdexcept>
#include <cmath>

using std::memcpy;
using std::ostream;
using std::invalid_argument;
using std::abs;

namespace mco {

Point::Point(unsigned int dimension) : dimension_(dimension) {
	values_ = new double[dimension];
	for(unsigned int i = 0; i < dimension; ++i)
		values_[i] = 0;
}

Point::Point(const double *values, unsigned int dimension)  : dimension_(dimension) {
	values_ = new double[dimension];
	memcpy(values_, values, dimension * sizeof(*values));
}

Point::Point(double value, unsigned int dimension) : dimension_(dimension) {
	values_ = new double[dimension];
	for(unsigned int i = 0; i < dimension; ++i)
		values_[i] = value;
}

Point::~Point() noexcept {
	delete[] values_;
}

Point::Point(const Point &p) : dimension_(p.dimension_) {
	values_ = new double[p.dimension_];
	memcpy(values_, p.values_, dimension_ * sizeof(*values_));
}

Point & Point::operator=(const Point &p) throw(invalid_argument) {
	if(dimension_ != p.dimension_)
		throw invalid_argument("Cannot assign to point with different dimension.");

	for(unsigned int i = 0; i < dimension_; ++i)
		values_[i] = p[i];

	return *this;
}

Point Point::operator+(const Point &p) const {
	double * values = new double[dimension_];
	for(unsigned int i = 0; i < dimension_; ++i)
		values[i] = (*this)[i] + p[i];

	Point point = Point(values, dimension_);
	delete[] values;
	return point;
}

Point Point::operator-() const {
	return (*this) * -1;
}


Point Point::operator-(const Point& p) const {
	return (*this) + -p;
}


Point & Point::operator+=(const Point &p) {
	for(unsigned int i = 0; i < dimension_; ++i)
		values_[i] += p[i];

	return *this;
}

Point & Point::operator-=(const Point &p) {
	for(unsigned int i = 0; i < dimension_; ++i)
		values_[i] -= p[i];

	return *this;
}

Point Point::operator*(double d) const {
	double * values = new double[dimension_];
	for(unsigned int i = 0; i < dimension_; ++i)
		values[i] = (*this)[i] * d;

	Point point = Point(values, dimension_);
	delete[] values;
	return point;
}

Point & Point::operator*=(double d) {
	for(unsigned int i = 0; i < dimension_; ++i)
		values_[i] *= d;

	return *this;
}

double Point::operator*(Point &point) const {
	double sum = 0;
	for(unsigned int i = 0; i < dimension_; ++i)
		sum += (*this)[i] * point[i];
	return sum;
}

bool Point::is_equal(const Point &p, double epsilon) const {
	for(unsigned int i = 0; i < dimension_; ++i)
		if(abs((*this)[i] - p[i]) > epsilon)
			return false;
	return true;
}

bool Point::is_not_equal(const Point &p, double epsilon) const {
	return !is_equal(p, epsilon);
}

bool Point::is_less_or_equal(const Point &p, double epsilon) const {
	for(unsigned int i = 0; i < dimension_; ++i)
		if((*this)[i] - p[i] > epsilon)
			return false;
	return true;
}

bool Point::is_less(const Point &p, double epsilon) const {
	for(unsigned int i = 0; i < dimension_; ++i)
		if((*this)[i] - p[i] > -epsilon)
			return false;
	return true;
}

bool Point::Dominates(const Point &p, double epsilon) const {
	bool equal = true;

	for(unsigned int i = 0; i < dimension_; i++) {
		if((*this)[i] - p[i] > epsilon)
			return false;

		if(abs((*this)[i] - p[i]) > epsilon)
			equal = false;
	}

	return !equal;
}

bool Point::is_lexicographic_less(const Point &point, double epsilon) const {
	for(unsigned int i = 0; i < dimension_; ++i) {
		if((*this)[i] - point[i] < -epsilon) {
			return true;
		} else if(point[i] - (*this)[i] < -epsilon) {
			return false;
		}
	}
	return false;
}

bool Point::is_lexicographic_less_or_equal(const Point &point, double epsilon) const {
	for(unsigned int i = 0; i < dimension_; ++i) {
		if((*this)[i] - point[i] < -epsilon) {
			return true;
		} else if((*this)[i] - point[i] > epsilon) {
			return false;
		}
	}
	return true;
}

double &Point::operator[](unsigned int index) {
	return values_[index];
}

double Point::operator[](unsigned int index) const {
	return values_[index];
}

unsigned int Point::dimension() const {
	return dimension_;
}

const double * Point::AsArray() const {
	return values_;
}

Point * Point::Null(unsigned int dimension) {
	return new Point(0.0, dimension);
}

Point * Point::One(unsigned int dimension) {
	return new Point(1.0, dimension);
}

std::ostream & operator<<(std::ostream &os, const Point &point) {
	os << "(";
	for(unsigned int i = 0; i < point.dimension_ - 1; ++i)
		os << point[i] << ", ";

	os << point[point.dimension_ - 1] << ")";

	return os;
}

bool LexicographicPointComparator::operator()(const Point * point1, const Point * point2) const {
	return point1->is_lexicographic_less(*point2, epsilon_);
}

}
