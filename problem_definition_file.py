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
data.mesh((0.0, 2.0, 0.0, 1.0, 0.0, 1.0), (65, 33, 33))

#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin = 1.0e-9
Emax = 1.0
penal = 3.0
data.material(Emin, Emax, 0.3, penal)

# filter: (type, radius)
# filter types: density = 1, sensitivity = 0
data.filter(1, 0.08)

# optimizer: (maxIter)
data.mma(4)

#volumefraction = 0.12
#nEl = data.nElements
#nnEl = nEl - (64 * 32)
#print(nEl)
#print((nEl - (64 * 32)))

#def roof(el):
 #   if el > nnel:
  #      val = 1.0
   # else:
   #     val = -1.0
   # return val

# passive elements
#data.passive(roof)

# BC: ('name' ,[checker], [setter])
bc1 = data.bc([1],[1])

# loadcases: (list of bc)
data.loadcases([bc1])


# Calculate the objective function
# objective input: (design variable value, SED)
def objective(xp, uKu):
    return (Emin + np.power(xp, penal) * (Emax - Emin)) * uKu

def sensitivity(xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

#def Constraint(xp, uKu):
 #   return (xp / nEl) - volumefraction

#def ConstraintSensitivity(xp, uKu):
 #   return 1.0 / nEl

# Callback implementation
data.obj(objective)
data.obj_sens(sensitivity)

# Volume constraint is standard, input (volume fraction)
#data.volumeConstraint(0.12)

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







# put in some list in a list: bc = [], [[bc 1],[bc 2],[bc 3],[bc 4],[bc 5],[bc 6]], 
# [ load case 1: [bc 1, bc 2], loadcase 2: [bc 3]]
# [ bc 1:[name, [checker], [setter]], [bc 2:...], [ bc etc...]], 
# checker:[lcoorp: i or i + 1, xc: 0 or 2] 
# setter:[vec:N or RHS, dof: i or i + 1, value:0.0 load in N] 
# parsing as a Python tuple

#data.bc(())
#data.bc((  [   ["fixed", [["lcoorp[i]", "xc[0]"]], [["N", "i", 0.0], ["N", "i + 1", 0.0], ["N", "i + 2", 0.0]],
#               [ "load", [["lcoorp[i]", "xc[1]"], ["lcoorp[i + 2]", "xc[4]"]], [["RHS", "i + 2", -0.001]]]
#            ]
#))
#bc1 = np.zeros((33*33))
#c = 0

# passing node numbers as numpy arrays
#for i in np.arange(0, data.nDOF, 1.):
#    if i == 0 or i % 100 == 0:
#        print(i, c)
#        bc1[c] = i
#        bc1[c+1] = i + 1
#        bc1[c+2] = i + 2 
#        c += 3  

#print(bc1)
#print(type(bc1))
#print(bc1.shape)
#print(bc1.dtype)
#print(bc1.nbytes)

#data.bc(bc1)