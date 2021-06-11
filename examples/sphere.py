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
# mesh: (domain: x, y, z, center)(mesh: number of nodes)
# 1/8
data.structuredGrid((0.0, 1.0, 0.0, 1.2, 0.0, 1.2), (161, 193, 193))

# readin STL file in binary format
# stl read: (encoding, backround, treshold, box around stl: (min corner)(max corner), full path to file)
# Passive elements: 1.0
# Active elements: -1.0
# Solid elements: 2.0
# Rigid elements: 3.0
# Void elements: 4.0
# Do not overwrite: 0.0
data.stlread(
    3.0,
    -1.0,
    4,
    (0.0, 0.0, 0.0),
    (1.0, 1.2, 1.2),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/rigidtorsion.stl",
)
data.stlread(
    4.0,
    0.0,
    4,
    (0.0, 0.0, 0.0),
    (1.0, 1.2, 1.2),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/torsionvoid.stl",
)

# Optional printing:
# print(data.nNodes)
# print(data.nElements)
# print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 1.0
data.material(Emin, Emax, nu, dens, penal)

# setup heavyside projection filter (betaFinal, betaInit, eta)
# data.projection(64.0, 1.0, 0.05)

# OR

# # robust formulation (betaFinal, betaInit, eta, delta)
data.robust(64.0, 1.0, 0.5, 0.1)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1,
data.filter(2, 0.065)

# optimizer: (maxIter, tol)
data.mma(2000, 0.01)

# loadcases: (# of loadcases)
data.loadcases(1)

# bc: (loadcase, type, [checker: dof index], [checker: values], [setter: dof index], [setter: values], parametrization)
# bc: (loadcase, type, [coordinate axis], [coordinate value], [coordinate axis], [bc value], parametrization)
# symmetry
data.bc(0, 1, [2], [0.0], [0, 1], [0.0, 0.0], 0)
data.bc(0, 1, [1], [0.0], [0, 2], [0.0, 0.0], 0)

# 1/8
data.bc(0, 1, [0], [1.0], [1, 2], [0.0, 0.0], 0)
data.bc(0, 2, [0, 1, 2], [0.0125, 0.0, 0.1625], [1], [0.01], 0)
data.bc(0, 2, [0, 1, 2], [0.0125, 0.1625, 0.0], [2], [-0.01], 0)

materialvolumefraction = 0.02
nEl = data.nael
rigidVol = data.nrel * 10.0


# Calculate the objective function
# objective input: (design variable value, SED)
def objective(comp, sumXp, volfrac):
    return comp


def sensitivity(xp, uKu, penal):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu


def constraint(comp, sumXp, volfrac):
    return (sumXp - rigidVol) / nEl - volfrac


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
