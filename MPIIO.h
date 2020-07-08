#ifndef MPIIO_H
#define MPIIO_H

// Include necessary libraries
#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include <mpi.h>
#include <petsc/private/dmdaimpl.h>
#include <petscdmda.h>
#include <string>
#include <vector>

/* -----------------------------------------------------------------------------
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
Copyright (C) 2013-2019,

This MPIIO implementation is licensed under Version 2.1 of the GNU
Lesser General Public License.

This MMA implementation is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This Module is distributed in the hope that it will be useful,implementation
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
-------------------------------------------------------------------------- */

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
    //void scale_shift(double scale_, V3 shift_);
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
	//GridBox(V3& minCorner_, V3& macCorner_, double dx_);
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


class MPIIO {
  public:
    // ------------- METHODS ------------------------------------------

    MPIIO(DM da_nodes, int nPfields, std::string pnames, int nCfields, std::string cnames);
    ~MPIIO();

    // NOT CLEAN INTERFACE: REPLACE BY STD::PAIR OR SUCH !!!!!!
    PetscErrorCode WriteVTK(DM da_nodes, Vec U, Vec x, Vec xTilde, Vec xPhys, PetscInt itr);

  private:
    // -------------- METHODS -----------------------------------------

    void abort(std::string errorMsg, std::string position);

    unsigned long int sum(unsigned long int* startPos, unsigned long int nel);
    // --------------- MEMBERS ----------------------------------------
    int         MPI_IS;               //!< The size of an MPI unsigned long integer in bytes
    int         MPI_FS;               //!< The size of an MPI float in bytes
    int         MPI_CS;               //!< The size of an MPI char in bytes
    int         nDom;                 //!< Number of domains
    int*        nPFields;             //!< Number of point fields in each domain
    int*        nCFields;             //!< Number of cell/element fields in each domain
    MPI_Offset  offset;               //!< The offset of each thread in the file
    int         rank;                 //!< The processor rank
    int         ncpu;                 //!< The number of cpus
    int         nodesPerElement;      //!< Number of nodes per element
    bool        firstFieldOutputDone; //!< Will be set to true after first field output
    std::string filename;             //!< Output filename
    MPI_File    fh;                   //!< Filehandle

    void Allocate(std::string info, const int nDom, const int nPFields[], const int nCFields[],
                  unsigned long int nPointsMyrank[], unsigned long int nCellsMyrank[],
                  unsigned long int nodesPerElement, std::string pFNames, std::string cFNames);
    // std::string filename = "/home/naage/PETSc/output2.dat");

    void writePoints(int domain, float coordinates[]);

    void writeCells(int domain, unsigned long int elements[], unsigned long int cellsOffset0[],
                    unsigned long int cellsTypes0[]);

    void writePointFields(unsigned long int timeStep, int domain, float fields[],
                          std::string newFilename = "notDefined");

    void writeCellFields(int domain, float fields[]);
    // ------------ MEMBERS  - maybe they can be private too -----------
    unsigned long int* nPoints;  //!< The number of points in each domain in each thread
    unsigned long int* nCells;   //!< The number of elements/cells in each domain in each thread
    unsigned long int* nPointsT; //!< The total number of points in each domain
    unsigned long int* nCellsT;  //!< The total number of cells in each domain

    // Converters needed for PETSc adaptation
    unsigned long int *nPointsMyrank, *nCellsMyrank;
    float *            workPointField, *workCellField;
    PetscErrorCode     DMDAGetElements_3D(DM dm, PetscInt* nel, PetscInt* nen, const PetscInt* e[]);
};
/** @example
        An illustrative example explaining how to use the class
        @code
  string userDefined = "something"
  int nDom 		 = 2;
  int nPFields[nDom]= {2, 1};
  string pFieldNames = "density, pdensity, density";
  int nCFields[nDom]= {1, 1};
  string cFieldNames = "density, density";
  int nPointsMyrank[nDom] = {10, 20}; // Arbitrary numbers
  int nCellsMyrank[nDom] = {9, 19};

  // Initialize
  MPIIO outputObject = DFEMMMPIIO(userDefined, nDom, nPFields, pFieldNames,
                                                                nCFields,
  cFieldNames, rank, nPoints, nCells);
  // Both the total number of points/cells and the corresponding
  // processor dependent numbers are determined in the constructor!
  // This does of course require MPI communication, but is only done once.


  // Output coordinates/points/vertices:
  // First domain:
  float points1[3*nPoints[0]];
  // ... put coordinates into the array points1
  outputObject.writePoints(0, points1);
  // Second domain:
  float points1[3*nPoints[1]];
  // ... put coordinates into the array points2
  outputObject.writePoints(1, points2);

  // Output elements:
  // First domain:
  int elements1[USER_SPECIFIED_SIZE];
  // ... put elements into the array elements1
  outputObject.writeCells(0, elements1);
  // Second domain
  int elements2[USER_SPECIFIED_SIZE];
  // ... put elements into the array
  outputObject.writeCells(1,  elements2);

  // The thought is that elements should be a long list with integers, where one
  integer
  // specifies the element type and the following integers specify point
  connectivity.

  // ...  Perform your computations

  // Output fields at a time point (timeStep)
  // First point fields are outputted, then cell fields

  // Output point field:
  // First domain:
  // Put all field variables into an array pFieldsD1, should of course
  // be in the same order as your point list, and the fields should come
  // after each other.
  outputObject.writePointField(timeStep, 0, pFieldsD1);
  // Second domain:
  // Put all field variables into an array pFieldsD2.
  outputObject.writePointField(1, pFieldsD2);

  // Output cell/element fields:
  // First domain:
  // Put all field variables into an array cFieldsD1, should of course
  // be in the same order as your element list, and the fields should come after
  each other. outputObject.writePointField(0, cFieldsD1);
  // Second domain:
  // Put all field variables into an array cFieldsD2.
  outputObject.writePointField(1, cFieldsD2);

  // Repeat as many times as you want.
        @endcode
*/

#endif
