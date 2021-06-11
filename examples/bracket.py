#!/usr/bin/env python3

# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

import numpy as np

import topoptlib

# step 1:
# Create data class to store input data
data = topoptlib.Data()

# step 2:
# define input data
# mesh: (domain: x, y, z, centers)(mesh: number of nodes)
data.structuredGrid((0.0, 192.0, 0.0, 64.0, 0.0, 104.0), (193, 65, 105))

# readin STL file in binary format
# stl read and voxelize: (encoding, backround, treshold, box around stl: (min corner)(max corner), full path to file)
# Passive elements: 1.0
# Active elements: -1.0
# Solid elements: 2.0
# Rigid elements: 3.0
# Void elements: 4.0
# Do not overwrite: 0.0
data.stlread(
    -1.0,
    1.0,
    8,
    (-23.0, -1.0, -103.0),
    (169.0, 63.0, 1.0),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineDesignDomainFine.stl",
)
data.stlread(
    2.0,
    0.0,
    4,
    (-23.0, -1.0, -103.0),
    (169.0, 63.0, 1.0),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineSolidDomainFine.stl",
)
data.stlread(
    3.0,
    0.0,
    8,
    (-23.0, -1.0, -103.0),
    (169.0, 63.0, 1.0),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineRigidDomainFine.stl",
)

# Optional printing:
# print(data.nNodes)
# print(data.nElements)
# print(data.nDOF)
# print(data.nael)
# print(data.nsel)
# print(data.nrel)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, Dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, Dens, penal)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1
data.filter(1, 5.0)

# optimizer: (maxIter, tol)
data.mma(400, 0.01)


def parametrization(lcx, lcy, lcz):
    # radius of bolt face d14.17= r7.5 mm
    # centers are
    # x 172.5, z 87.5
    # x 17.5, z 49.5
    # x 166.5, z 35.5
    # x 22.5, z 87.5
    val = 1.0
    if np.power((lcx - 172.5) / 7.5, 2) + np.power((lcz - 87.5) / 7.5, 2) < 1:
        val = 0.0
        # print('1', lcx, lcz)

    if np.power((lcx - 17.5) / 7.5, 2) + np.power((lcz - 49.5) / 7.5, 2) < 1:
        val = 0.0
        # print('2', lcx, lcz)

    if np.power((lcx - 166.5) / 7.5, 2) + np.power((lcz - 35.5) / 7.5, 2) < 1:
        val = 0.0
        # print('3', lcx, lcz)

    if np.power((lcx - 22.5) / 7.5, 2) + np.power((lcz - 87.5) / 7.5, 2) < 1:
        val = 0.0
        # print('4', lcx, lcz)

    return val


# loadcases: (# of loadcases)
data.loadcases(4)

# use parametrization fuction
data.bcpara(parametrization)

# bc: (loadcase, type, [checker: dof index], [checker: values], [setter: dof index], [setter: values], parametrization)
# bc: (loadcase, type, [coordinate axis], [coordinate value], [coordinate axis], [bc value], parametrization)
# up
data.bc(0, 1, [1], [1.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(0, 1, [1], [7.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(0, 2, [0, 1, 2], [89.0, 45.0, 21.0], [1], [0.008], 0)

# out
data.bc(1, 1, [1], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(1, 1, [1], [7.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(1, 2, [0, 1, 2], [89.0, 45.0, 21.0], [2], [-0.0085], 0)

# 42deg
data.bc(2, 1, [1], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(2, 1, [1], [7.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(2, 2, [0, 1, 2], [89.0, 45.0, 21.0], [2], [-0.0095 * np.sin(np.deg2rad(42))], 0)
data.bc(2, 2, [0, 1, 2], [89.0, 45.0, 21.0], [1], [0.0095 * np.cos(np.deg2rad(42))], 0)

# torsion
data.bc(3, 1, [1], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(3, 1, [1], [7.0], [0, 1, 2], [0.0, 0.0, 0.0], 1)
data.bc(3, 2, [0, 1, 2], [68.0, 45.0, 23.0], [2], [0.005], 0)
data.bc(3, 2, [0, 1, 2], [111.0, 45.0, 19.0], [2], [-0.005], 0)

materialvolumefraction = 0.3
nEl = data.nael
rigidVol = data.nrel * 10.0
solidVol = data.nsel * 1.0


# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp, volfrac):
    return comp


def sensitivity(xp, uKu, penal):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu


def constraint(comp, sumXp, volfrac):
    return (sumXp - rigidVol - solidVol) / nEl - materialvolumefraction


def constraintSensitivity(xp, uKu, penal):
    return 1.0 / nEl


# Callback implementation
data.obj(objective)
data.objsens(sensitivity)

# Define constraint
data.cons(constraint)
data.conssens(constraintSensitivity)

# Homogeniuos initial condition
data.initialcondition(materialvolumefraction)

# Output vtr files
data.vtr(20)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()
