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
data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0), (129, 65, 65))

# Optional printing:
# print(data.nNodes)
# print(data.nElements)
# print(data.nDOF)

# material: (Emin, Emax, nu, penal)
Emin, Emax, nu, Dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, Dens, penal)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1
# using 0.08, 0.04 or 0.02
data.filter(1, 0.08)

# optimizer: (maxIter, tol)
data.mma(4000, 0.01)

# loadcases: (# of loadcases)
data.loadcases(2)

# bc: (loadcase, type, [checker: dof index], [checker: values], [setter: dof index], [setter: values], parametrization)
# bc: (loadcase, type, [coordinate axis], [coordinate value], [coordinate axis], [bc value], parametrization)
data.bc(0, 1, [0], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 2, [0, 2], [2.0, 0.0], [2], [-0.001], 0)
data.bc(0, 2, [0, 1, 2], [2.0, 0.0, 0.0], [2], [-0.0005], 0)
data.bc(0, 2, [0, 1, 2], [2.0, 1.0, 0.0], [2], [-0.0005], 0)

data.bc(1, 1, [0], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(1, 2, [0, 1], [2.0, 1.0], [1], [0.001], 0)
data.bc(1, 2, [0, 1, 2], [2.0, 1.0, 0.0], [1], [0.0005], 0)
data.bc(1, 2, [0, 1, 2], [2.0, 1.0, 1.0], [1], [0.0005], 0)

materialvolumefraction = 0.24
nEl = data.nElements


# Calculate the objective function, senitivity, constraint and constraint sensitivity
def objective(comp, sumXp, volfrac):
    return comp


def sensitivity(xp, uKu, penal):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu


def constraint(comp, sumXp, volfrac):
    return sumXp / nEl - materialvolumefraction


def constraintSensitivity(xp, uKu, penal):
    return 1.0 / nEl


# Callback implementatio
data.obj(objective)
data.objsens(sensitivity)

# Define constraint
data.cons(constraint)
data.conssens(constraintSensitivity)

# Local volume constraint input: (Rlocvol, alpha)
data.localVolume(0.16, 0.12)

# Homogeniuos initial condition
data.initialcondition(materialvolumefraction)

# Output vtr files
data.vtr(20)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()
