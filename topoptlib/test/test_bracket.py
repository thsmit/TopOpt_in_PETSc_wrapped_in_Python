#!/usr/bin/python3

# Author: Thijs Smit, May 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

##EULER_MEMORY="2500"
#NCPU=32
#WALL_TIME="24:00"

import topoptlib
import numpy as np

# step 1:
# Create data class to store input data
data = topoptlib.Data()

# step 2:
# define input data
# mesh: (domain: x, y, z, center)(mesh: number of nodes)
data.structuredGrid((0.0, 192.0, 0.0, 64.0, 0.0, 104.0, 1.0, 0.0, 0.0, 0.0, 0.0), (193, 65, 105))
#data.structuredGrid((0.0, 192.0, 0.0, 64.0, 0.0, 104.0, 1.0, 0.0, 0.0), (385, 129, 209))

# readin STL file in binary format
# stl read: ((box around stl: (min corner)(max corner))full path to file)
data.stlread_domain((-23.0, -1.0, -103.0), (169.0, 63.0, 1.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineDesignDomainFine.stl')

# stl read: load the solid domain into the same coordinate system as the design domain
data.stlread_solid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineSolidDomainFine.stl')

# stl read: load the solid domain into the same coordinate system as the design domain
data.stlread_rigid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineRigidDomainFine.stl')

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, Dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, Dens, penal)

# Iterpro = 20, stepsize = 0.125, penalinit = 1.0, penalfinal = 3.0
#data.continuation()
#data.projection()

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, 
data.filter(1, 5.0)

# optimizer: (maxIter)
data.mma(11)

def parametrization(lcx, lcy, lcz):
    # radius of bolt face 14.17= r7.5 mm
    # centers are 
    # x 172.5, z 87.5
    # x 17.5, z 49.5
    # x 166.5, z 35.5
    # x 22.5, z 87.5
    val = 1.0
    if np.power( (lcx - 172.5) / 7.5, 2) + np.power( (lcz - 87.5) / 7.5, 2) < 1:
        val = 0.0
        #print('1', lcx, lcz)

    if np.power( (lcx - 17.5) / 7.5, 2) + np.power( (lcz - 49.5) / 7.5, 2) < 1:
        val = 0.0
        #print('2', lcx, lcz)

    if np.power( (lcx - 166.5) / 7.5, 2) + np.power( (lcz - 35.5) / 7.5, 2) < 1:
        val = 0.0
        #print('3', lcx, lcz)

    if np.power( (lcx - 22.5) / 7.5, 2) + np.power( (lcz - 87.5) / 7.5, 2) < 1:
        val = 0.0
        #print('4', lcx, lcz)

    return val

# loadcases: (# of loadcases)
data.loadcases(4)

# use parametrization fuction
data.bcpara(parametrization)

# bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values], parametrization)
# up
data.bc(0, 1, [1, 6], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(0, 2, [0, 89.0, 1, 45.0, 2, 21.0], [1], [0.008], 2)

# out
data.bc(1, 1, [1, 6], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(1, 2, [0, 89.0, 1, 45.0, 2, 21.0], [2], [-0.0085], 2)

#42deg
data.bc(2, 1, [1, 6], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(2, 2, [0, 89.0, 1, 45.0, 2, 21.0], [2], [-0.0095*np.sin(np.deg2rad(42))], 2)
data.bc(2, 2, [0, 89.0, 1, 45.0, 2, 21.0], [1], [0.0095*np.cos(np.deg2rad(42))], 2)

# torsion
data.bc(3, 1, [1, 6], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(3, 2, [0, 68.0, 1, 45.0, 2, 23.0], [2], [0.005], 2)
data.bc(3, 2, [0, 111.0, 1, 45.0, 2, 19.0], [2], [-0.005], 2)

materialvolumefraction = 0.3
#nEl = data.nElements
nEl = data.nael - data.nrel - data.nsel
#nEl = data.nael
rigidVol = data.nrel * 10.0
solidVol = data.nsel * 1.0

# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp, xp, uKu):
    return comp

def sensitivity(sumXp, xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

def constraint(comp, sumXp, xp, uKu):
    #print('sumXp',sumXp)
    #print('rigidVol',rigidVol)
    #print('solidVol',solidVol)
    #print('nEl',nEl)
    return (sumXp - rigidVol - solidVol) / nEl - materialvolumefraction

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
# post processing, generate .vtu file to be viewed in paraview
#if complete:
#    data.vtu()

