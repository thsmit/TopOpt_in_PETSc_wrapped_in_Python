#include "V3.h"
#include <cmath>

V3::V3() {
}

V3::V3(char* facet) {
	char f1[4] = { facet[0], facet[1], facet[2], facet[3] };
	char f2[4] = { facet[4], facet[5], facet[6], facet[7] };
	char f3[4] = { facet[8], facet[9], facet[10], facet[11] };
	float xx = *((float*)f1);
	float yy = *((float*)f2);
	float zz = *((float*)f3);
	x = double(xx);
	y = double(yy);
	z = double(zz);
}

V3::V3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {
}

V3::~V3() {
}

ostream& operator<<(ostream& os, V3 v) {
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}

const V3 V3::operator+(const V3& rv) const {
	return V3(x + rv.x, y + rv.y, z + rv.z);
}

const V3 V3::operator-(const V3& rv) const {
	return V3(x - rv.x, y - rv.y, z - rv.z);
}

const V3 V3::operator*(const double& rv) const {
	return V3(x * rv, y * rv, z * rv);
}

const V3 V3::operator/(const double& rv) const {
	return V3(x / rv, y / rv, z / rv);
}

const V3 V3::operator^(const V3& rv) const {
	return V3(
		y * rv.z - z * rv.y,
		z * rv.x - x * rv.z,
		x * rv.y - y * rv.x);
}

const double V3::operator&(const V3& rv) const {
	return x * rv.x + y * rv.y + z * rv.z;
}

V3& V3::operator+=(const V3& rv) {
	x += rv.x;
	y += rv.y;
	z += rv.z;
	return *this;
}

V3& V3::operator-=(const V3& rv) {
	x -= rv.x;
	y -= rv.y;
	z -= rv.z;
	return *this;
}

V3& V3::operator*=(const double& rv) {
	x *= rv;
	y *= rv;
	z *= rv;
	return *this;
}

V3& V3::operator/=(const double& rv) {
	x /= rv;
	y /= rv;
	z /= rv;
	return *this;
}

double V3::normal() const {
	return sqrt(x * x + y * y + z * z);
}

V3 V3::zero = V3(0.0, 0.0, 0.0);
V3 V3::one = V3(1.0, 1.0, 1.0);
