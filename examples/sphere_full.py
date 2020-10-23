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

# Note:
# set nlvls to 1 in TopOpt.cc and LinearElasticity.cc

# step 2:
# define input data
# mesh: (domain: x, y, z, center)(mesh: number of nodes)
#data.structuredGrid((0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.5, 1.0, 1.5, 0.0, 0.0), (5, 5, 5)) # full
#data.structuredGrid((0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0), (3, 3, 3)) # 1/8
data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.5, 1.0, 0.0, 0.0), (5, 3, 3)) # 1/4

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, dens, penal = 1.0e-9, 1.0, 0.3, 1.0, 1.0
data.material(Emin, Emax, nu, dens, penal)

#data.continuation()
#data.projection()

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(1, 0.000025)

# optimizer: (maxIter)
data.mma(1) # or 1

# loadcases: (# of loadcases)
data.loadcases(1)

# bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values])
#data.bcpara(parametrization)

# full
'''
data.bc(0, 1, [0, 1, 1, 6, 2, 6], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 6, 2, 7], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 6, 2, 8], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 7, 2, 6], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 7, 2, 7], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 7, 2, 8], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 8, 2, 6], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 8, 2, 7], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 8, 2, 8], [0, 1, 2], [0.0, 0.0, 0.0], 0)

data.bc(0, 2, [0, 0, 1, 8, 2, 8], [1], [0.1*np.sin(np.deg2rad(45))], 0)
data.bc(0, 2, [0, 0, 1, 8, 2, 8], [2], [-0.1*np.sin(np.deg2rad(45))], 0)

data.bc(0, 2, [0, 0, 1, 6, 2, 8], [1], [0.1*np.sin(np.deg2rad(45))], 0)
data.bc(0, 2, [0, 0, 1, 6, 2, 8], [2], [0.1*np.sin(np.deg2rad(45))], 0)

data.bc(0, 2, [0, 0, 1, 6, 2, 6], [1], [-0.1*np.sin(np.deg2rad(45))], 0)
data.bc(0, 2, [0, 0, 1, 6, 2, 6], [2], [0.1*np.sin(np.deg2rad(45))], 0)

data.bc(0, 2, [0, 0, 1, 8, 2, 6], [1], [-0.1*np.sin(np.deg2rad(45))], 0)
data.bc(0, 2, [0, 0, 1, 8, 2, 6], [2], [-0.1*np.sin(np.deg2rad(45))], 0)
'''


# 1/8

data.bc(0, 1, [0, 1, 1, 0, 2, 7], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 7, 2, 7], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [0, 1, 1, 7, 2, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 1, [2, 0], [0, 1], [0.0, 0.0], 0)
data.bc(0, 1, [1, 0], [0, 2], [0.0, 0.0], 0)
data.bc(0, 2, [0, 0, 1, 7, 2, 7], [1], [0.1*np.sin(np.deg2rad(45))], 0)
data.bc(0, 2, [0, 0, 1, 7, 2, 7], [2], [-0.1*np.sin(np.deg2rad(45))], 0)



materialvolumefraction = 1.0
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

# Homogeniuos initial condition
data.initialcondition(materialvolumefraction)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()