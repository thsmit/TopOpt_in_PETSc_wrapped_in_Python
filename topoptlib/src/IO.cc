// Import necessary stuff
#include "IO.h"
#include <cstdlib> // To get the exit function


// STL reader
// reference: https://github.com/zhulianhua/Voxelizer
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>

// CGAL
//#include <CGAL/Simple_cartesian.h>

// VTK
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkCellArray.h>
//#include <vtkCellData.h>
//#include <vtkDataArray.h>
//#include <vtkDataSetMapper.h>
//#include <vtkIdList.h>
//#include <vtkNamedColors.h>
//#include <vtkPointData.h>
//#include <vtkPoints.h>
//#include <vtkPolyhedron.h>
//#include <vtkProperty.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkXMLUnstructuredGridWriter.h>

using namespace std;

void writeVTKPolyData() {
	
	// reading binary file
	//printf("Reading .dat\n");
	std::ifstream FIN("output_00000.dat", std::ios::in | std::ios::binary);
    
	if (FIN)
	{
		// print mesh
		std::string header;
		std::getline(FIN, header);
		std::cout << header << std::endl;

		// load info
		// # Dom, # points, # cells, #pFields, #cFields, #nodes per element
		//int nDom, nPointsT, nCellsT, nPFields, nCFields, nodesPerElement;
		
		//FIN.read(nDom, sizeof(nDom));

		//float a, b, c, d, e, f;

    	//FIN >> nDom;

    	//printf("nDom: %i\t", nDom);

		//cout << "Integer: " << nDom << endl;

		// point field names
		//std::string pointFieldNames;
		//std::getline(FIN, pointFieldNames);
		//std::cout << pointFieldNames << std::endl;

		// cell field names
		//std::string cellFieldNames;
		//std::getline(FIN, cellFieldNames);
		//std::cout << cellFieldNames << std::endl;

		

		//while (FIN) {
		//	std::string strInput;
		//	std::getline(FIN, strInput);
		//	std::cout << strInput << std::endl;
		//}
		//FIN.read(headInfo, 26);
		//std::cout << "INFO: " << headInfo << std::endl;
		//char nTriRaw[4];
		//stlFile.read(nTriRaw, 4);
		//unsigned numTri = *((unsigned*)nTriRaw);
		//triangles.resize(numTri);
		//for (auto&& tri : triangles)
		//{
		//	char triRaw[50];
		//	stlFile.read(triRaw, 50);
		//	V3 norm(triRaw);
		//	V3 p1(triRaw+12);
		//	V3 p2(triRaw+24);
		//	V3 p3(triRaw+36);
		//	tri = Triangle(norm, p1, p2, p3);
		//}
	}
	else
	{
		std::cerr << ".data file openning error!" << std::endl;
	}

	// writing VTKdata
	
	// create polyhedron (cube)
	//vtkIdType pointIds[8] = {0, 1, 2, 3, 4, 5, 6, 7};

	//vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	//points->InsertNextPoint(-1.0,-1.0,-1.0);
	//points->InsertNextPoint( 1.0,-1.0,-1.0);
	//points->InsertNextPoint( 1.0, 1.0,-1.0);
	//points->InsertNextPoint(-1.0, 1.0,-1.0);
	//points->InsertNextPoint(-1.0,-1.0, 1.0);
	//points->InsertNextPoint( 1.0,-1.0, 1.0);
	//points->InsertNextPoint( 1.0, 1.0, 1.0);
	//points->InsertNextPoint(-1.0, 1.0, 1.0);

	//vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	//vtkIdType face0[4] = {0, 3, 2, 1};
	//vtkIdType face1[4] = {0, 4, 7, 3};
	//vtkIdType face2[4] = {4, 5, 6, 7};
	//vtkIdType face3[4] = {5, 1, 2, 6};
	//vtkIdType face4[4] = {0, 1, 5, 4};
	//vtkIdType face5[4] = {2, 3, 7, 6};

	//faces->InsertNextCell(4, face0);
	//faces->InsertNextCell(4, face1);
	//faces->InsertNextCell(4, face2);
	//faces->InsertNextCell(4, face3);
	//faces->InsertNextCell(4, face4);
	//faces->InsertNextCell(4, face5);

	//vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//ugrid->SetPoints(points);
	//ugrid->InsertNextCell(VTK_POLYHEDRON, 8, pointIds, 6, faces->GetPointer());

	// Here we write out the cube.
	//vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	//writer->SetInputData(ugrid);
	//writer->SetFileName("polyhedron.vtu");
	//writer->SetDataModeToAscii();
	//writer->Update();
	
	//printf("Close .dat, but keep the file as backup\n");
}

// reference: www.jgxsoft.com
bool checkBinarySTL(char * buffer) {

	// Determine if a file is binary STL.
	bool bbinary = true;
	size_t spnsz, spnsz0;

	// Look for the first non-space character.
	spnsz = strspn(buffer, " ");

	char ctstr[6];  // Enough space to hold "solid\0" and "facet\0".
	// Copy the first five characters from the location of the first non-space
	// character to ctstr.
	//strncpy_s(ctstr, &buffer[spnsz], 5);
	strncpy(ctstr, &buffer[spnsz], 5);

	ctstr[5] = '\0';
	char csolid[] = "solid\0";

	// If this is an ASCII STL file, then the first string should be "solid".
	if (!strcmp(ctstr, csolid)) {
		// This file might be binary or text. To be certain we need to do a further test.
		// Read past the next new line. If this is a text file, there should be a newline.
		// The next token should be 'facet'.

		spnsz0 = 5 + spnsz;

		char * pch = strchr(&buffer[spnsz0], '\n');  // Look for the first instance of '\n'.
		// If pch is NULL then this is a binary STL file.
		if (pch) {
			pch++;

			spnsz = strspn(pch, " "); // Look for the first instance not of ' '.
			spnsz0 = spnsz;

			spnsz = strcspn(pch + spnsz0, " "); // Look for the first instance of ' '.

			if (spnsz == 5) {
				// Check for 'facet'.
				//strncpy_s(ctstr, pch + spnsz0, 5);
				strncpy(ctstr, pch + spnsz0, 5);
				ctstr[5] = '\0';

				char cfacet[] = "facet\0";
				if (!strcmp(ctstr, cfacet)) {
					// This file is beyond reasonable doubt ASCII STL.
					bbinary = false;
				}
			}
		}
	}

	return(bbinary);


}

//#define DEBUG
V3::V3() { }

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

V3::V3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) { }

V3::~V3() { }

std::ostream& operator<<(std::ostream& os, V3 v) {
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

Triangle::Triangle() { }

Triangle::Triangle(V3 p1_, V3 p2_, V3 p3_) : p1(p1_), p2(p2_), p3(p3_)
{
	//TODO:
	// Compute norm via p1, p2, p3
}

Triangle::Triangle(V3 norm_, V3 p1_, V3 p2_, V3 p3_) :
	norm(norm_), p1(p1_), p2(p2_), p3(p3_)
{}

Triangle::~Triangle() { }

std::ostream& operator<<(std::ostream& os, Triangle& tri)
{
	os << "(\n  " << tri.p1 << "\n  " << tri.p2 <<"\n  " << tri.p3 << "\n)";
	return os;
}

// global functions
// determin z plane / triangle relative position type
std::vector<int> tri_plane(Triangle& tri_, double z_){
	std::vector<int> ans;
	if (tri_.p1.z > z_) ans.push_back(1);
	if (tri_.p2.z > z_) ans.push_back(2);
	if (tri_.p3.z > z_) ans.push_back(3);
	return ans;
}

// get the section line of triangle with a z plane
Line tri_sec_plane(Triangle& tri_, std::vector<int>& aboveP_, double z_) {
	Line line1, line2;
	if (aboveP_.size() == 1) {
		if (aboveP_[0] == 1) {
			line1 = Line(tri_.p1, tri_.p2);
			line2 = Line(tri_.p1, tri_.p3);
		}
		else if (aboveP_[0] == 2) {
			line1 = Line(tri_.p2, tri_.p1);
			line2 = Line(tri_.p2, tri_.p3);
		}
		else {
			line1 = Line(tri_.p3, tri_.p1);
			line2 = Line(tri_.p3, tri_.p2);
		}
	}
	else {  // 2 above points
		// find the below point first
		if (aboveP_[0] == 1 && aboveP_[1] == 2) {
			line1 = Line(tri_.p3, tri_.p1);
			line2 = Line(tri_.p3, tri_.p2);
		} 
		else if (aboveP_[0] == 1 && aboveP_[1] == 3) {
			line1 = Line(tri_.p2, tri_.p3);
			line2 = Line(tri_.p2, tri_.p1);
		}
		else {
			line1 = Line(tri_.p1, tri_.p2);
			line2 = Line(tri_.p1, tri_.p3);
		}
	}
	//
	V3 p1, p2;
	p1 = line1.p_cross_z_plane(z_);
	p2 = line2.p_cross_z_plane(z_);
	return Line(p1, p2);
}

std::ostream& operator<<(std::ostream& os_, Bbox& box_)
{
	os_ << "Bbox(" << std::endl 
		<< "minCorner: \t" << box_.minCorner << std::endl
		<< "maxCorner: \t" << box_.maxCorner << std::endl
		<< "Extend: \t" << box_.get_extend() << std::endl << ")" << std::endl;
	return os_;
}

Line::Line(V3& p1_, V3& p2_) : p1(p1_), p2(p2_) {};

V3 Line::p_cross_z_plane(double z) const {
	V3 ans = p1 + (p2 - p1) * (z - p1.z) / (p2.z - p1.z + 1e-10);
	return ans;
}

V3 Line::p_cross_y_plane(double y) const {
	V3 ans = p1 + (p2 - p1) * (y - p1.y) / (p2.y - p1.y + 1e-30);
	return ans;
}

V3 Line::p_cross_x_plane(double x) const {
	V3 ans = p1 + (p2 - p1) * (x - p1.x) / (p2.x - p1.x);
	return ans;
}

Geometry::Geometry(std::string fname) {
	//std::cout << "Reading geometry ..." << std::endl;
	read_stl_file(fname);
	//std::cout << "Num of triangle " << get_num_tri() << std::endl;
	set_bound();
}

Geometry::~Geometry() {}

GridBox::GridBox(V3& minCorner_, V3& maxCorner_, double dx_) 
	: Bbox(minCorner_, maxCorner_), dx(dx_)
{
	// the size in x/y/z direction has to be multiple of dx
	//std::cout << "Generate grid ..." << std::endl;
	V3 extend = maxCorner - minCorner;
	int nxI = (int)(extend.x / dx);
	int nyI = (int)(extend.y / dx);
	int nzI = (int)(extend.z / dx);
	double xRes = extend.x - nxI * dx;
	double yRes = extend.y - nyI * dx;
	double zRes = extend.z - nzI * dx;
	if (xRes > 1e-6 * dx || yRes > 1e-6 * dx || zRes > 1e-6 * dx)
	{
		//std::cerr << "Grid box size is not dx*Nxyz" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	gridNum = int3(nxI, nyI, nzI);
}

GridBox::GridBox(V3& minCorner_, double dx_, int3 gridNum_)
	: dx(dx_), gridNum(gridNum_)
{
	double lx = dx_ * gridNum.nx;
	//std::cout << "lx " << lx << std::endl;
	double ly = dx_ * gridNum.ny;
	double lz = dx_ * gridNum.nz;
	minCorner = minCorner_;
	maxCorner = minCorner + V3(lx, ly, lz);
}

void Geometry::read_stl_file(std::string fname) {
	std::ifstream stlFile(fname, std::ios::in | std::ios::binary);
    
    char headInfo[80] = "";
	// read 80 byte header
	if (stlFile)
	{
		stlFile.read(headInfo, 80);
		//std::cout << "STL file comment " << headInfo << std::endl;
		char nTriRaw[4];
		stlFile.read(nTriRaw, 4);
		unsigned numTri = *((unsigned*)nTriRaw);
		triangles.resize(numTri);
		for (auto&& tri : triangles)
		{
			char triRaw[50];
			stlFile.read(triRaw, 50);
			V3 norm(triRaw);
			V3 p1(triRaw+12);
			V3 p2(triRaw+24);
			V3 p3(triRaw+36);
			tri = Triangle(norm, p1, p2, p3);
		}
	}
	else
	{
		std::cerr << "File openning error!" << std::endl;
	}
}

void Geometry::scale_shift(double scale_, V3 shift_) {
	// keep a track of the current shift and scale v.s. to the STL file
	for (auto&& tri : triangles)
	{
		// scale then shift
		tri.p1 *= scale_;
		tri.p1 += shift_;
		tri.p2 *= scale_;
		tri.p2 += shift_;
		tri.p3 *= scale_;
		tri.p3 += shift_;
		tri.norm *= scale_;
	}
	// min/max corner of the bound
	V3 minCorner = bound.get_minCorner();
	V3 maxCorner = bound.get_maxCorner();
	minCorner *= scale_;
	minCorner += shift_;
	maxCorner *= scale_;
	maxCorner += shift_;
	bound = Bbox(minCorner, maxCorner);
}

void Geometry::set_bound() {
	double minX, minY, minZ, maxX, maxY, maxZ;
	minX = minY = minZ =  1.0e30;
	maxX = maxY = maxZ = -1.0e30;
	for (auto&& tri : triangles)
	{
		minX = (tri.p1.x < minX) ? tri.p1.x : minX;
		minX = (tri.p2.x < minX) ? tri.p2.x : minX;
		minX = (tri.p3.x < minX) ? tri.p3.x : minX;
		minY = (tri.p1.y < minY) ? tri.p1.y : minY;
		minY = (tri.p2.y < minY) ? tri.p2.y : minY;
		minY = (tri.p3.y < minY) ? tri.p3.y : minY;
		minZ = (tri.p1.z < minZ) ? tri.p1.z : minZ;
		minZ = (tri.p2.z < minZ) ? tri.p2.z : minZ;
		minZ = (tri.p3.z < minZ) ? tri.p3.z : minZ;

		maxX = (tri.p1.x > maxX) ? tri.p1.x : maxX;
		maxX = (tri.p2.x > maxX) ? tri.p2.x : maxX;
		maxX = (tri.p3.x > maxX) ? tri.p3.x : maxX;
		maxY = (tri.p1.y > maxY) ? tri.p1.y : maxY;
		maxY = (tri.p2.y > maxY) ? tri.p2.y : maxY;
		maxY = (tri.p3.y > maxY) ? tri.p3.y : maxY;
		maxZ = (tri.p1.z > maxZ) ? tri.p1.z : maxZ;
		maxZ = (tri.p2.z > maxZ) ? tri.p2.z : maxZ;
		maxZ = (tri.p3.z > maxZ) ? tri.p3.z : maxZ;
	}
	V3 minCorner = V3(minX, minY, minZ);
	V3 maxCorner = V3(maxX, maxY, maxZ);
	bound = Bbox(minCorner, maxCorner);
}

Voxelizer::Voxelizer(Geometry& geo_, GridBox& grid_) : geo(geo_), grid(grid_) {
	// scale (in the z direction) and shift the geometry to fit the grid
	//std::cout << "Generating voxilzer ..." << std::endl;
	Bbox bound = geo.get_bound();
	//std::cout << "GEO bound = " << bound << std::endl;

	//double dx = grid.get_dx();
	//std::cout << "dx vox= " << dx << std::endl;
	//double scale = (grid.get_extend().z - 2.0 * dx )/ bound.get_extend().z;
	// double scale = 1.126; // small MBB Example
	//double scale = 1.05; // big MBB Example
	//std::cout << "scale vox= " << scale << std::endl;
	//geo.scale_shift(scale, V3::zero);


	//V3 shift = grid.get_minCorner() - geo.get_bound().get_minCorner();
	//shift += V3(1.0, 0.5, 0.5); // MBB Example
	//std::cout << "shift vox= " << shift << std::endl;
	//geo.scale_shift(1.0, shift);

	int3 gridNum = grid.get_gridNum();
	// + 1 added
	int nx = gridNum.nx + 1;
	int ny = gridNum.ny + 1;
	int nz = gridNum.nz + 1;
	//std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
	long numTotal = nx * ny * nz;
	
	// flag is point data!
	flag = new char[numTotal];
	for (int i = 0; i < numTotal; i++)
		flag[i] = 0;

	for (int iz = 1; iz <= nz; iz++)
	{
		std::vector<Line> lines;
		get_z_sections(lines, iz);
		for (int iy = 1; iy <= ny; iy++)
		{
			std::vector<int> xids;
			get_xid_cross(xids, iy, lines, iz);
			int n_change = 0;
			bool isBlack = false;
			if (xids.empty()) {
				for (int ix = 0; ix < nx; ix++) {
					flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
				}
			} 
			else{
				for (int ix = 0; ix < nx; ix++) {
					//cout << "ix = " << ix << "   n_change  = " << n_change << endl;
					if (n_change < xids.size() && ix == xids[n_change] ) {
						n_change++;
						if (isBlack) {
							isBlack = !isBlack;
							flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
						}
						else {
							flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
							isBlack = !isBlack;
						}
					}
					else{
						flag[(iz - 1) * nx * ny + (iy - 1) * nx + ix] = isBlack;
					}
				}
			}
		}
	}
}

Voxelizer::~Voxelizer()
{
	//std::cout << "Deleting the flag memory!" << std::endl;
	delete []flag;
}

void Voxelizer::get_relevant_triangles(std::vector<Triangle>& tri_, int iz_) const {
	// case 0: 3 points below, 0 points above
	// case 1: 2 points below, 1 points above
	// case 2: 1 points below, 2 points above
	// case 3: 0 points below, 3 points above
	for (int i = 0; i < geo.get_num_tri(); i++) {
		Triangle triTmp = geo.get_tri(i);
		double dx = grid.get_dx();
		double zmin = grid.get_minCorner().z;
		std::vector<int> aboveP = tri_plane(triTmp, iz_ * dx + zmin);
		if (aboveP.size() > 0 && aboveP.size() < 3) {
			tri_.push_back(triTmp);
		}
	}
}

void Voxelizer::get_z_sections(std::vector<Line>& lines_, int iz_, std::vector<Triangle>& tris_) const {}

void Voxelizer::get_z_sections(std::vector<Line>& lines_, int iz_) const {
	lines_.clear();
	double z = grid.get_dx() * iz_ + grid.get_minCorner().z;
	for (int i = 0; i < geo.get_num_tri(); i++) {
		Triangle triTmp = geo.get_tri(i);
		std::vector<int> aboveP = tri_plane(triTmp, z);
		if (aboveP.size() > 0 && aboveP.size() < 3) {
			Line lineTmp = tri_sec_plane(triTmp, aboveP, z);
			lines_.push_back(lineTmp);
			std::cout.precision(16);
#ifdef DEBUG
			if (lines_.size() == 56 && iz_ == 4) {
				cout << "Abnominal line end : " << endl;
				cout << lineTmp.p2 << endl;
				cout << "Above P size " << aboveP.size() << endl;
				cout << "Plane Z  " << z << endl;
				cout << "Tri id is " << i << endl;
				cout << "Tri = " << endl;
				cout << triTmp << endl;
			}
#endif
		}
	}
}

void Voxelizer::get_xid_cross(std::vector<int>& xids_, int iy_, std::vector<Line>& lines_, int iz_) const {
	V3 minCorner = grid.get_minCorner();
	double dx = grid.get_dx();
	double y = iy_ * dx + minCorner.y;
	xids_.clear();
	for (auto && lineTmp : lines_) 
		if ((lineTmp.p1.y - y) * (lineTmp.p2.y - y) < 0.0) { // line cross y
			V3 cross = lineTmp.p_cross_y_plane(y);
			int xx = int((cross.x - minCorner.x) / dx);
			xids_.push_back(xx);
		}
	// sort
	std::sort(xids_.begin(), xids_.end());
}

void Voxelizer::write_vtk_image() {
	//std::cout << "Writing vtk image ... " << std::endl;
	int3 gridNum = grid.get_gridNum();
	int nx = gridNum.nx;
	int ny = gridNum.ny;
	int nz = gridNum.nz;

	std::ofstream of("flag.vtk");
	of << "# vtk DataFile Version 2.0" << std::endl;
	of << "Flag data" << std::endl;
	of << "ASCII" << std::endl << std::endl;
	of << "DATASET STRUCTURED_POINTS" << std::endl;
	of << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
	of << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
	of << "SPACING " << 1 << " " << 1 << " " << 1 << std::endl << std::endl;
	of << "POINT_DATA " << nx * ny * nz << std::endl;
	of << "SCALARS flag int" << std::endl;
	of << "LOOKUP_TABLE deafult" << std::endl;
	int n = nx * ny * nz;
	for (int i = 0; i < n; i++)
		of << (int)flag[i] << std::endl;
}

