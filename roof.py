#!/usr/bin/python3

# Author: Thijs Smit, April 2020
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
# mesh: (domain)(mesh: number of nodes)
dx = 3.0
nx = 73
nz = 25
data.mesh((0.0, dx, 0.0, 0.5, 0.0, 1.0), (nx, 13, nz))

# material: (Emin, Emax, nu, penal)
Emin = 1.0e-9
Emax = 1.0
penal = 3.0
data.material(Emin, Emax, 0.3, penal)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, denity + heaviside projection = 2
rmin = (dx/float(nx))*2
print(rmin)
data.filter(3, rmin)

# projection:

# optimizer: (maxIter)
data.mma(500)

nEl = data.nElements
nnEl = nEl - ((nx - 1) * (nz - 1))

def roof(el):
    if el > nnel:
        val = 1.0
    else:
        val = -1.0
    return val

# passive elements
data.passive(roof)

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(xp, uKu):
    return (Emin + np.power(xp, penal) * (Emax - Emin)) * uKu

def sensitivity(xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

# Callback implementation
data.obj(objective)
data.obj_sens(sensitivity)

# Volume constraint is standard, input (volume fraction)
data.volumeConstraint(0.25)

# Define additional constraints
data.const(Constraint)
data.const_sens(ConstraintSensitivity)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()

#print(data.it)
#print(data.trueFx)
#print(data.scaledFx)

# step 4:
# post processing, generate .vtu file to be viewed in paraview
if complete:
    data.vtu()

# step 5:
# generate .stl file from final design for 3D printing
#if complete:
#    data.stl()