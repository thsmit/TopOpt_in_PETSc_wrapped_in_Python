#ifndef IO_H
#define IO_H

/*
Author: Thijs Smit, Dec 2020
Copyright (C) 2020 ETH Zurich

Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

#include <string>
#include <vector>

void writeVTKPolyData();
bool checkBinarySTL(char *);

// STL reader
class V3
{
public:
	V3();
	V3(char* facet);
	V3(double x_, double y_, double z_);
	~V3();
	friend std::ostream& operator<<(std::ostream& os, V3 v);
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

class Triangle
{
public:
	Triangle();
	Triangle(V3 p1_, V3 p2_, V3 p3_);
	Triangle(V3 norm_, V3 p1_, V3 p2_, V3 p3_);
	~Triangle();
	friend std::ostream& operator<<(std::ostream& os, Triangle& tri);
	V3 p1, p2, p3, norm;
};

class int3 {
public:
	int3() : nx(0), ny(0), nz(0) {};
	int3(int nx_, int ny_, int nz_) : nx(nx_), ny(ny_), nz(nz_) {};
	int nx, ny, nz;
};

class Line {
public:
	Line() {};
	Line(V3& p1_, V3& p2_);
	V3 get_length() const;
	V3 p_cross_y_plane(double y) const;
	V3 p_cross_z_plane(double z) const;
	V3 p_cross_x_plane(double x) const;
	V3 p1;
	V3 p2;
};

class Bbox {
public:
	Bbox() {};
	Bbox(V3& minCorner_, V3& maxCorner_) : minCorner(minCorner_), maxCorner(maxCorner_) {};
	inline V3 get_minCorner() const { return minCorner; };
	inline V3 get_maxCorner() const { return maxCorner; };
	inline V3 get_extend() const { return V3(maxCorner - minCorner); };
	friend std::ostream& operator<<(std::ostream& os_, Bbox& box_);
	~Bbox() {};
protected:
	V3 minCorner;
	V3 maxCorner;
};

class Geometry {
  public:
	  Geometry(std::string fname_);
	  inline int get_num_tri() const { return (int)triangles.size(); };
	  //inline V3 get_refPoint() const { return refPoint; };
    inline Bbox get_bound() const { return bound; };
    void set_bound();
    //void set_refPoint(V3& refPoint_);
    void scale_shift(double scale_, V3 shift_);
    inline Triangle get_tri(int i) const { return triangles[i];};
    ~Geometry();
  private:
    void read_stl_file(std::string fname);
    V3 refPoint;
    std::vector<Triangle> triangles;
    Bbox bound;
};

class GridBox : public Bbox {
public:
	GridBox(V3& minCorner_, V3& macCorner_, double dx_);
	GridBox(V3& minCorner_, double dx_, int3 gridNum_);
	inline double get_dx() const { return dx; };
	inline int3 get_gridNum() const { return gridNum; };
private:
	double dx;
	int3 gridNum;
};

class Voxelizer {
public:	
	// constructor do the vexelization
	Voxelizer(Geometry& geo_, GridBox& grid_);
	// return the voxelized flag;
	const char* get_flag() const { return flag; };
	void write_vtk_image();
	// ~descruction
	~Voxelizer();
private:
	// get triangles that intersect with plane iz_
	void get_relevant_triangles(std::vector<Triangle>& tri_,  int iz_) const;

	// TO BE DELETED
	// get intersection of the geometry with the plane iz_;
	void get_z_sections(std::vector<Line>& lines_, int iz_, std::vector<Triangle>& tris_) const;

	void get_z_sections(std::vector<Line>& lines_, int iz_) const;

	// get xid of the cells that intersection with lines at iy_
	void get_xid_cross(std::vector<int>& xids_, int iy_, std::vector<Line>& lines_) const;
	// get xid of the cells that intersection with lines at iy_  for debug
	void get_xid_cross(std::vector<int>& xids_, int iy_, std::vector<Line>& lines_, int iz_) const;
	// flag data
	char* flag;


	// geometry
	Geometry geo;
	// grid box
	GridBox grid;
};

#endif