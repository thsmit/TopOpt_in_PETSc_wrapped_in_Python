#ifndef V3_H
#define V3_H
#include <iostream>
using namespace std;
class V3
{
public:
	V3();
	V3(char* facet);
	V3(double x_, double y_, double z_);
	~V3();
	friend ostream& operator<<(ostream& os, V3 v);
	const V3 operator+(const V3& rv) const;
	const V3 operator-(const V3& rv) const;
	const V3 operator*(const double& rv) const;
	const V3 operator/(const double& rv) const;
	// extern product
	const V3 operator^(const V3& rv) const;
	// inner product
	const double operator&(const V3& rv) const;
	// add
	V3& operator+=(const V3& rv);
	// substract
	V3& operator-=(const V3& rv);
	// multiply
	V3& operator*=(const double& rv);
	// divide
	V3& operator/=(const double& rv);
	// static variables
	static V3 zero;
	static V3 one;
	// normal vector
	double normal() const;
	// coordinate
	double x, y, z;
};
#endif
