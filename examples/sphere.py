#!/usr/bin/env python3

# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

import topoptlib
import numpy as np

# step 1:
# Create data class to store input data
data = topoptlib.Data()

# step 2:
# define input data
# mesh: (domain: x, y, z, center)(mesh: number of nodes)
# 1/8
data.structuredGrid((0.0, 1.0, 0.0, 1.2, 0.0, 1.2, 0.1625, 0.0125, 0.1875, 0.9875, 0.0), (161, 193, 193))

# 1/4
#data.structuredGrid((0.0, 2.0, 0.0, 1.2, 0.0, 1.2, 0.1625, 0.0125, 0.1875, 0.9875, 0.0), (161, 97, 97))


# readin STL file in binary format
# TO DO: allow for ASCII format
# stl read: ((box around stl: (min corner)(max corner)), full path to file)
# 1/8
data.stlread(3.0, -1.0, 4, (0.0, 0.0, 0.0), (1.0, 1.2, 1.2), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/rigidtorsion.stl')
data.stlread(4.0, 0.0, 4, (0.0, 0.0, 0.0), (1.0, 1.2, 1.2), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/torsionsolid.stl')

# 1/4
#data.stlread(3.0, -1.0, 4, (0.0, 0.0, 0.0), (2.0, 1.2, 1.2), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/torsionrigid14.stl')
#data.stlread(4.0, 0.0, 4, (0.0, 0.0, 0.0), (2.0, 1.2, 1.2), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/torsionvoid14.stl')

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 1.0
data.material(Emin, Emax, nu, dens, penal)

# setup continuation of penalization: (Pinitial, Pfinal, stepsize, IterProg)
#data.continuation(1.0, 3.0, 1.0, 100)

# setup heavyside projection filter (betaFinal, stepsize, eta) update of beta every 20 iterations
data.projection(64.0, 1.0, 0.5) 

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(1, 0.025)

# optimizer: (maxIter, tol)
data.mma(1500, 0.001)

# loadcases: (# of loadcases)
data.loadcases(1)

# symmetry
data.bc(0, 1, [2, 0], [0, 1], [0.0, 0.0], 0)
data.bc(0, 1, [1, 0], [0, 2], [0.0, 0.0], 0)

# 1/8
data.bc(0, 1, [0, 1], [1, 2], [0.0, 0.0], 0)
data.bc(0, 2, [0, 0, 1, 0, 2, 6], [1], [0.01], 0)
data.bc(0, 2, [0, 0, 1, 6, 2, 0], [2], [-0.01], 0)

# back # 1/4
#data.bc(0, 1, [0, 1, 1, 0, 2, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
#data.bc(0, 1, [0, 1, 1, 0, 2, 6], [0, 1, 2], [0.0, 0.0, 0.0], 0)
#data.bc(0, 1, [0, 1, 1, 6, 2, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)

materialvolumefraction = 0.02
nEl = data.nael
rigidVol = data.nrel * 10.0

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp):
    return comp

def sensitivity(xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

def constraint(comp, sumXp):
    return (sumXp - rigidVol) / nEl - materialvolumefraction

def constraintSensitivity(xp, uKu):
    return 1.0 / nEl

# Callback implementation
data.obj(objective)
data.objsens(sensitivity)

# Define constraint
data.cons(constraint)
data.conssens(constraintSensitivity)

# Homogeniuos initial condition
data.initialcondition(materialvolumefraction)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()

# step 4:
# generate .vtu file to be viewed in paraview
# generate .x3d file for import and modification in Blender
#if complete:
#    data.vtu()
#    data.x3d()

# step 5:
# post processing by smoothning and output binary .stl
# if complete:
#    smooth = data.smoothening()
# if smooth:
#    data.stl(smooth)
#    data.stp(smooth)

