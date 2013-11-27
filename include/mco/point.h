#pragma once
/*
 * point.h
 *
 *  Created on: 15.03.2013
 *      Author: fritz
 */

#ifndef POINT_H_
#define POINT_H_

#include <ostream>
#include <stdexcept>

namespace mco {

class LexicographicPointComparator;

class Point {

	double * values_ = nullptr;
	const unsigned int dimension_ = 0;

public:

	Point() = delete;
	~Point() noexcept;
	explicit Point(unsigned int dimension);
	Point(const double *values, unsigned int dimension);
	Point(double value, unsigned int dimension);

	Point(const Point &p);
	Point & operator=(const Point &p) throw(std::invalid_argument);

	Point operator+(const Point &p) const;
	Point operator-(const Point &p) const;
	Point operator-() const;
	Point & operator+=(const Point &p);
	Point & operator-=(const Point &p);
	Point operator*(double d) const;
	Point & operator*=(double d);
	double operator*(Point &point) const;

	bool is_equal(const Point &p, double epsilon) const;
	bool is_not_equal(const Point &p, double epsilon) const;
	bool is_less_or_equal(const Point &p, double epsilon) const;
	bool is_less(const Point &p, double epsilon) const;

	bool is_lexicographic_less(const Point &point, double epsilon) const;
	bool is_lexicographic_less_or_equal(const Point &point, double epsilon) const;

	bool operator==(const Point &d) = delete;
	bool operator!=(const Point &d) = delete;
	bool operator<=(const Point &d) = delete;
	bool operator>=(const Point &d) = delete;
	bool operator<(const Point &d) = delete;
	bool operator>(const Point &d) = delete;

	double &operator[](unsigned int index);
	double operator[](unsigned int index) const;
	unsigned int dimension() const;

	bool Dominates(const Point &p, double epsilon) const;

	const double * AsArray() const;

	static Point * Null(const unsigned int dimension);
	static Point * One(const unsigned int dimension);

	friend std::ostream & operator<<(std::ostream &, const Point &);
	friend LexicographicPointComparator;

};

std::ostream & operator<<(std::ostream &os, const Point &point);

class LexicographicPointComparator {
public:
	LexicographicPointComparator() = delete;

	LexicographicPointComparator(double epsilon) : epsilon_(epsilon) {}

	bool operator()(const Point * point1, const Point * point2) const;

private:
	double epsilon_;
};

}


#endif /* POINT_H_ */
