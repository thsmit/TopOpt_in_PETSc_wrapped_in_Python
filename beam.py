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
#data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (17, 9, 9))
#data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (33, 17, 17))
data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (65, 33, 33))

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, Dens, penal = 1.0e-9, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, Dens, penal)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(2, 0.08)

# optimizer: (maxIter)
data.mma(40)

# loadcases: (# of loadcases)
data.loadcases(1)

# bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values])
data.bc(0, 1, [0, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 2, [0, 1, 2, 4], [2], [-0.001], 0)
data.bc(0, 2, [0, 1, 1, 2, 2, 4], [2], [-0.0005], 0)
data.bc(0, 2, [0, 1, 1, 3, 2, 4], [2], [-0.0005], 0)

#data.bc(1, 1, [[0, 0]], [0, 1, 2], [0.0, 0.0, 0.0], None)
#data.bc(1, 2, [[0, 1], [1, 3]], [1], [0.001], None)
#data.bc(1, 2, [[0, 1], [1, 3], [2, 4]], [1], [0.0005], None)
#data.bc(1, 2, [[0, 1], [1, 3], [2, 5]], [1], [0.0005], None)

materialvolumefraction = 0.12
nEl = data.nElements

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp, xp, uKu):
    return comp

def sensitivity(sumXp, xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

def constraint(comp, sumXp, xp, uKu):
    return sumXp / nEl - materialvolumefraction

def constraintSensitivity(sumXp, xp, uKu):
    return 1.0 / nEl

# Callback implementation
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