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
#data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (129, 65, 65))
#data.structuredGrid((0.0, 4.0, 0.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0), (17, 9, 9))
#data.structuredGrid((0.0, 4.0, 0.0, 2.0, 0.0, 2.0, 1.125, 3.125, 0.375, 0.5, 1.5), (33, 17, 17))
data.structuredGrid((0.0, 4.0, 0.0, 2.0, 0.0, 2.0, 1.03125, 3.03125, 0.46875, 0.5, 1.5), (129, 65, 65))
#data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (65, 33, 33))

data.stlread_domain((0.0, 0.0, 0.0), (4.0, 2.0, 2.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/beam.stl')
data.stlread_solid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/beamrigid.stl')
data.stlread_rigid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/input/beamrigid.stl')

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, Dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, Dens, penal)

# setup continuation of penalization: (Pinitial, Pfinal, stepsize)
#data.continuation()

# setup heavyside projection filter (betaFinal, stepsize, eta)
#data.projection() 

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(1, 0.08)

# optimizer: (maxIter)
data.mma(1600)

# loadcases: (# of loadcases)
data.loadcases(1)

# bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values])
data.bc(0, 1, [0, 6], [0, 1, 2], [0.0, 0.0, 0.0], 0)
#data.bc(0, 2, [0, 7, 2, 8], [2], [-0.001], 0)
#data.bc(0, 1, [0, 7, 1, 9, 2, 8], [1], [0.0], 0)
#data.bc(0, 1, [0, 7, 1, 10, 2, 8], [1], [0.0], 0)
data.bc(0, 2, [0, 7, 1, 9, 2, 8], [2], [-0.0005], 0)
data.bc(0, 2, [0, 7, 1, 10, 2, 8], [2], [-0.0005], 0)

#data.bc(1, 1, [[0, 0]], [0, 1, 2], [0.0, 0.0, 0.0], None)
#data.bc(1, 2, [[0, 1], [1, 3]], [1], [0.001], None)
#data.bc(1, 2, [[0, 1], [1, 3], [2, 4]], [1], [0.0005], None)
#data.bc(1, 2, [[0, 1], [1, 3], [2, 5]], [1], [0.0005], None)

materialvolumefraction = 0.12
nEl = data.nael
solidVol = data.nsel * 1.0
rigidVol = data.nrel * 10.0
print('nEl', data.nael)

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp, xp, uKu):
    return comp

def sensitivity(sumXp, xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

def constraint(comp, sumXp, xp, uKu):
    return (sumXp - solidVol - rigidVol) / nEl - materialvolumefraction

def constraintSensitivity(sumXp, xp, uKu):
    return 1.0 / nEl

# Callback implementatio
data.obj(objective)
data.objsens(sensitivity)

# Define constraint
data.cons(constraint)
data.conssens(constraintSensitivity)

# Volume constraint is standard, input (volume fraction)
data.initialcondition(materialvolumefraction)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()

# step 4:
# post processing, generate .vtu file to be viewed in paraview
#if complete:
#    data.vtu()