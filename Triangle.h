#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "V3.h"
class Triangle
{
public:
	Triangle();
	Triangle(V3 p1_, V3 p2_, V3 p3_);
	Triangle(V3 norm_, V3 p1_, V3 p2_, V3 p3_);
	~Triangle();
	friend ostream& operator<<(ostream& os, Triangle& tri);
	V3 p1, p2, p3, norm;
};

#endif