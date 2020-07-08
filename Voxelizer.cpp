#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "V3.h"
#include "Voxelizer.h"

//#define DEBUG

using namespace std;

// global functions
// determin z plane / triangle relative position type
vector<int> tri_plane(Triangle& tri_, double z_){
	vector<int> ans;
	if (tri_.p1.z > z_) ans.push_back(1);
	if (tri_.p2.z > z_) ans.push_back(2);
	if (tri_.p3.z > z_) ans.push_back(3);
	return ans;
}

// get the section line of triangle with a z plane
Line tri_sec_plane(Triangle& tri_, vector<int>& aboveP_, double z_) {
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

ostream& operator<<(ostream& os_, Bbox& box_)
{
	os_ << "Bbox(" << endl 
		<< "minCorner: \t" << box_.minCorner << endl
		<< "maxCorner: \t" << box_.maxCorner << endl
		<< "Extend: \t" << box_.get_extend() << endl << ")" << endl;
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

Geometry::Geometry(string fname) {
	cout << "Reading geometry ..." << endl;
	read_stl_file(fname);
	cout << "Num of triangle " << get_num_tri() << endl;
	set_bound();
}

Geometry::~Geometry() {}

GridBox::GridBox(V3& minCorner_, V3& maxCorner_, double dx_) 
	: Bbox(minCorner_, maxCorner_), dx(dx_)
{
	// the size in x/y/z direction has to be multiple of dx
	cout << "Generate grid ..." << endl;
	V3 extend = maxCorner - minCorner;
	int nxI = (int)(extend.x / dx);
	int nyI = (int)(extend.y / dx);
	int nzI = (int)(extend.z / dx);
	double xRes = extend.x - nxI * dx;
	double yRes = extend.y - nyI * dx;
	double zRes = extend.z - nzI * dx;
	if (xRes > 1e-6 * dx || yRes > 1e-6 * dx || zRes > 1e-6 * dx)
	{
		cerr << "Grid box size is not dx*Nxyz" << endl;
		exit(EXIT_FAILURE);
	}
	gridNum = int3(nxI, nyI, nzI);
}

GridBox::GridBox(V3& minCorner_, double dx_, int3 gridNum_)
	: dx(dx_), gridNum(gridNum_)
{
	double lx = dx_ * gridNum.nx;
	double ly = dx_ * gridNum.ny;
	double lz = dx_ * gridNum.nz;
	minCorner = minCorner_;
	maxCorner = minCorner + V3(lx, ly, lz);
}

void Geometry::read_stl_file(string fname) {
	ifstream stlFile(fname.c_str(), ios::in | ios::binary);
	char headInfo[80] = "";
	// read 80 byte header
	if (stlFile)
	{
		stlFile.read(headInfo, 80);
		cout << "STL file comment " << headInfo << endl;
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
		cerr << "File openning error!" << endl;
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
	cout << "Generating voxilzer ..." << endl;
	Bbox bound = geo.get_bound();
	cout << "GEO bound = " << bound << endl;

	double dx = grid.get_dx();
	//double scale = (grid.get_extend().z - 2.0 * dx )/ bound.get_extend().z;
	//geo.scale_shift(scale, V3::zero);

	//V3 shift = grid.get_minCorner() - geo.get_bound().get_minCorner();
	//shift += V3(2*dx, 2*dx, 2*dx);
	//geo.scale_shift(1.0, shift);

	int3 gridNum = grid.get_gridNum();
	int nx = gridNum.nx;
	int ny = gridNum.ny;
	int nz = gridNum.nz;
	cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << endl;
	long numTotal = nx * ny * nz;
	flag = new char[numTotal];
	for (int i = 0; i < numTotal; i++)
		flag[i] = 0;

	for (int iz = 1; iz <= nz; iz++)
	{
		vector<Line> lines;
		get_z_sections(lines, iz);
		for (int iy = 1; iy <= ny; iy++)
		{
			vector<int> xids;
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
	cout << "Deleting the flag memory!" << endl;
	delete []flag;
}

void Voxelizer::get_relevant_triangles(vector<Triangle>& tri_, int iz_) const {
	// case 0: 3 points below, 0 points above
	// case 1: 2 points below, 1 points above
	// case 2: 1 points below, 2 points above
	// case 3: 0 points below, 3 points above
	for (int i = 0; i < geo.get_num_tri(); i++) {
		Triangle triTmp = geo.get_tri(i);
		double dx = grid.get_dx();
		double zmin = grid.get_minCorner().z;
		vector<int> aboveP = tri_plane(triTmp, iz_ * dx + zmin);
		if (aboveP.size() > 0 && aboveP.size() < 3) {
			tri_.push_back(triTmp);
		}
	}
}

void Voxelizer::get_z_sections(vector<Line>& lines_, int iz_, vector<Triangle>& tris_) const {}

void Voxelizer::get_z_sections(vector<Line>& lines_, int iz_) const {
	lines_.clear();
	double z = grid.get_dx() * iz_ + grid.get_minCorner().z;
	for (int i = 0; i < geo.get_num_tri(); i++) {
		Triangle triTmp = geo.get_tri(i);
		vector<int> aboveP = tri_plane(triTmp, z);
		if (aboveP.size() > 0 && aboveP.size() < 3) {
			Line lineTmp = tri_sec_plane(triTmp, aboveP, z);
			lines_.push_back(lineTmp);
			cout.precision(16);
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

void Voxelizer::get_xid_cross(vector<int>& xids_, int iy_, vector<Line>& lines_, int iz_) const {
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
	cout << "Writing vtk image ..." << endl;
	int3 gridNum = grid.get_gridNum();
	int nx = gridNum.nx;
	int ny = gridNum.ny;
	int nz = gridNum.nz;
	ofstream of("flag.vtk");
	of << "# vtk DataFile Version 2.0" << endl;
	of << "Flag data" << endl;
	of << "ASCII" << endl << endl;
	of << "DATASET STRUCTURED_POINTS" << endl;
	of << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
	of << "ORIGIN " << 0 << " " << 0 << " " << 0 << endl;
	of << "SPACING " << 1 << " " << 1 << " " << 1 << endl << endl;
	of << "POINT_DATA " << nx * ny * nz << endl;
	of << "SCALARS flag int" << endl;
	of << "LOOKUP_TABLE deafult" << endl;
	int n = nx * ny * nz;
	for (int i = 0; i < n; i++)
		of << (int)flag[i] << endl;
}
