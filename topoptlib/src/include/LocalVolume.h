#ifndef __LOCALVOLUME__
#define __LOCALVOLUME__

#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include <petsc/private/dmdaimpl.h>
#include <iostream>
#include <math.h>

#include <topoptlib.h>

/* -----------------------------------------------------------------------------
Authors: Niels Aage, June 2019
Copyright (C) 2013-2019,

This LocalVolume implementation is licensed under Version 2.1 of the GNU
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


class LocalVolume{
  
public:
    // Constructor
    LocalVolume(DM da_nodes, Vec x, DataObj data); 
    
    // Destructor
    ~LocalVolume(); 
    
    // LocalVolume design variables
    PetscErrorCode Constraint(Vec x, PetscScalar *gx, Vec dx);
    
    // Get field  
    Vec GetField(){return xVol;};   
    
private:
  
    // Standard density/sensitivity filter matrix
    PetscScalar R;      // Radius for averaging
    PetscScalar pnorm;  // expnontent for p-norm
    PetscScalar alpha;  // max average value
    
    Vec xVol;
    Mat H; 		    // LocalVolume matrix
    Vec Hs; 		// LocalVolume "sum weight" (normalization factor) vector   
    
    // Mesh used for standard filtering
    DM da_elem;  	// da for image-filter field mesh
    
    // Setup datastructures for the filter
    PetscErrorCode SetUp(DM da_nodes);
    
    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);
    
};

#endif
