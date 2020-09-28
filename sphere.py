#!/usr/bin/python3

# Author: Thijs Smit, May 2020
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
#data.structuredGrid((0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0625, 1.9375, 0.0), (65, 65, 65))
#data.structuredGrid((0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0625, 1.9375, 0.0), (65, 65, 65))
#data.structuredGrid((0.0, 1.0, 0.0, 1.2, 0.0, 1.2, 0.2, 0.00625, 0.20625), (161, 193, 193))
#data.structuredGrid((0.0, 1.0, 0.0, 1.2, 0.0, 1.2, 0.2, 0.025, 0.2125), (81, 97, 97))
data.structuredGrid((0.0, 1.0, 0.0, 1.2, 0.0, 1.2, 0.2, 0.025, 0.2125), (161, 193, 193))


# readin STL file in binary format
# TO DO: allow for ASCII format
# stl read: ((box around stl: (min corner)(max corner)), full path to file)
data.stlread_domain((0.0, 0.0, 0.0), (1.0, 1.2, 1.2), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/domaintorsion.stl')
#data.stlread_domain((0.0, 0.0, 0.0), (1.0, 1.2, 1.2), 'full')
#data.stlread_domain((0.0, 0.0, 0.0), (1.2, 1.2, 1.2), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/hollowsphere.stl')
#data.stlread_domain((0.0, 0.0, 0.0), (1.0, 1.0, 1.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/tritext_demo.stl')
# data.stlread_domain('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/bracket/jetEngineDesignDomainFine.stl')
#data.stlread_solid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/solidtorsion.stl')
data.stlread_rigid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/rigidtorsionbig.stl')

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 1.0
data.material(Emin, Emax, nu, dens, penal)

#data.continuation()
#data.projection()

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(1, 0.025)

# optimizer: (maxIter)
data.mma(200) # or 1

# loadcases: (# of loadcases)
data.loadcases(1)

# define BC parametrization
# input: coordinates of node
# output: boolean, True or False
#def parametrization(lcx, lcy, lcz):
#    cenx = 1.0
#    ceny = 1.0
#    cenz = 1.0
#    xx = lcx - cenx
#    yy = lcy - ceny
#    zz = lcz - cenz

#    if np.power((yy / 0.1), (2)) + np.power((zz / 0.1), (2)) < 1:
#        if lcx == 0.0625:
#            print(lcx, lcy, lcz)

#    return np.power( np.power((yy / ceny), (2)) + np.power((zz / cenz), (2)), (1 / 0.10))

# bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values])
#data.bcpara(parametrization)
#data.bc(0, 1, [0, 6], [0, 1, 2], [0.0, 0.0, 0.0], 1)
#data.bc(0, 5, [0, 7], [1, 2], [-0.001], 1)
data.bc(0, 1, [0, 1], [0, 1, 2], [0.0, 0.0, 0.0], 0)
#data.bc(0, 1, [2, 4], [2], [0.0], 0)
#data.bc(0, 1, [1, 2], [1], [0.0], 0)
data.bc(0, 2, [0, 7, 1, 7, 2, 6], [1], [-0.001], 0)
data.bc(0, 2, [0, 7, 2, 7, 1, 6], [2], [0.001], 0)

#data.bc(1, 1, [[0, 0]], [0, 1, 2], [0.0, 0.0, 0.0], None)
#data.bc(1, 2, [[0, 1], [1, 3]], [1], [0.001], None)
#data.bc(1, 2, [[0, 1], [1, 3], [2, 4]], [1], [0.0005], None)
#data.bc(1, 2, [[0, 1], [1, 3], [2, 5]], [1], [0.0005], None)

materialvolumefraction = 0.02
#nEl = data.nElements
nEl = data.nael
rigidVol = data.nrel * 1000.0

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp, xp, uKu):
    return comp

def sensitivity(sumXp, xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

def constraint(comp, sumXp, xp, uKu):
    return (sumXp - rigidVol) / nEl - materialvolumefraction

def constraintSensitivity(sumXp, xp, uKu):
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

