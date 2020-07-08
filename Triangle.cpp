#include "Triangle.h"

Triangle::Triangle()
{
}

Triangle::Triangle(V3 p1_, V3 p2_, V3 p3_) : p1(p1_), p2(p2_), p3(p3_)
{
	//TODO:
	// Compute norm via p1, p2, p3
}

Triangle::Triangle(V3 norm_, V3 p1_, V3 p2_, V3 p3_) :
	norm(norm_), p1(p1_), p2(p2_), p3(p3_)
{}

Triangle::~Triangle()
{
}

ostream& operator<<(ostream& os, Triangle& tri)
{
	os << "(\n  " << tri.p1 << "\n  " << tri.p2 <<"\n  " << tri.p3 << "\n)";
	return os;
}
