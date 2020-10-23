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
# mesh: (domain: x, y, z)
#       (mesh: number of nodes: x, y, z)
#data.structuredGrid((0.0, 48.0, 0.0, 40.0, 0.0, 16.0), (97, 81, 33))
data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (17, 9, 9))

#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin = 1.0e-9
Emax = 1.0
penal = 3.0
data.material(Emin, Emax, 0.3, penal)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(0, 0.2)

# optimizer: (maxIter)
data.mma(2)

volumefraction = 0.12
nEl = data.nElements
#nnEl = nEl - (64 * 32)
#print(nEl)
#print((nEl - (64 * 32)))

def roof(el):
    #cx = 48
    #cy = 23
    #r = 5

    #x = lcx - cx
    #y = lcy - cy

    #print("el", el)

    #return (np.power((x / r), 2) + np.power((y / cy), 2))
    if el < (nEl - (64 * 32)):
        return -1.0
    else:
        return -1.0


# passive elements (parametrization function)
# only available for sensitivity and density filter, not for PDE filter
#data.passive(roof)

# loadcases: (# of loadcases)
data.loadcases(1)

# bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values])
data.bc(0, 1, [0, 0], [0, 1, 2], [0.0, 0.0, 0.0], None)
data.bc(0, 2, [0, 1, 2, 4], [2], [-0.001], None)
data.bc(0, 2, [0, 1, 1, 2, 2, 4], [2], [-0.0005], None)
data.bc(0, 2, [0, 1, 1, 3, 2, 4], [2], [-0.0005], None)

#data.bc(1, 1, [0, 0], [0, 1, 2], [0.0, 0.0, 0.0], None)
#data.bc(1, 2, [0, 1, 1, 3], [1], [0.001], None)
#data.bc(1, 2, [0, 1, 1, 3, 2, 4], [1], [0.0005], None)
#data.bc(1, 2, [0, 1, 1, 3, 2, 5], [1], [0.0005], None)

def param(lcx, lcy, lcz):
    cenx = (48.0 / 2) + 1
    ceny = (40.0 / 2) + 1
    cenz = (16.0 / 2) + 1
    xx = lcx - cenx
    yy = lcy - ceny
    zz = lcz - cenz

    return np.power( np.power((xx / cenx), (2 / 1)) + np.power((yy / ceny), (2 / 1)), (1 / 0.1))

# bc_para: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values], param)
#data.bc(2, 1, [2, 4], [0, 1, 2], [0.0, 0.0, 0.0], None)
#data.bc(2, 2, [2, 5], [2], [-0.001], param)

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(xp, uKu): 
    return (Emin + np.power(xp, penal) * (Emax - Emin)) * uKu

def sensitivity(xp, uKu): # make more general, check with Niels
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

# Callback implementation
data.obj(objective)
data.obj_sens(sensitivity)

# Volume constraint is standard, input (volume fraction)
data.volumeConstraint(volumefraction)

#def Constraint(xp, uKu):
 #   return (xp / nEl) - volumefraction

#def ConstraintSensitivity(xp, uKu):
 #   return 1.0 / nEl

# Define additional constraints
#data.const(Constraint)
#data.const_sens(ConstraintSensitivity)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()

#print(data.it)
#print(data.trueFx)
#print(data.scaledFx)

# step 4:
# post processing, generate .vtu file to be viewed in paraview
#if complete:
#    data.vtu()

# step 5:
# generate .stl file from final design for 3D printing
#if complete:
#    data.stl()