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
# mesh: (domain: x, y, z)(mesh: number of nodes)
data.structuredGrid((0.0, 84.0, 0.0, 70.0, 0.0, 112.0), (241, 201, 321))

# readin STL file in binary format
# stl read: (encoding, backround, threshold, box around stl: (min corner)(max corner), full path to file)
# Passive elements: 1.0
# Active elements: -1.0
# Solid elements: 2.0
# Rigid elements: 3.0
# Do not overwrite: 0.0
data.stlread(
    2.0,
    1.0,
    4,
    (16.0, 85.0, 1585.0),
    (100.0, 155.0, 1697.0),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/femur.stl",
)

# readin STL file in binary format
# stl read: (encoding, backround, threshold, box around stl: (min corner)(max corner), full path to file)
# Passive elements: 1.0
# Active elements: -1.0
# Solid elements: 2.0
# Rigid elements: 3.0
# Do not overwrite: 0.0
data.stlread(
    -1.0,
    0.0,
    8,
    (16.0, 85.0, 1585.0),
    (100.0, 155.0, 1697.0),
    "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/femur.stl",
)

# Optional printing:
# print(data.nNodes) # number of nodes
# print(data.nElements) # total number of elements, full domain
# print(data.nDOF) # total number of dofs, full domain
# print(data.nael) # number of active elements, design domain
# print(data.nsel) # number of solid elements, passive
# print(data.nrel) # number of rigid elements, passive

# material: (Emin, Emax, nu, dens, penal)
Emin, Emax, nu, dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, dens, penal)

# setup heavyside projection filter (betaFinal, betaInit, eta)
data.projection(64.0, 1.0, 0.5)

# filter: (type, radius)
# filter types: sensitivity = 0, density = 1, pde = 2,
data.filter(1, 0.5)

# optimizer: (maxIter, tol)
data.mma(2, 0.01)

# loadcases: (# of loadcases)
data.loadcases(1)

# bc: (loadcase, type, [coordinate axis], [coordinate value], [coordinate axis], [bc value], parametrization)
data.bc(
    0, 1, [2], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 0
)  # fully constraint at mid-diaphysis
data.bc(
    0,
    2,
    [0, 1, 2],
    [21.125, 22.75, 104.0],
    [0, 1, 2],
    [-2317 * np.sin(np.deg2rad(24)), 1.0, -2317 * np.cos(np.deg2rad(24))],
    0,
)  # Hip-joint contact force with 30deg
data.bc(
    0,
    2,
    [0, 1, 2],
    [50.375, 53.625, 84.5],
    [0, 1, 2],
    [-703 * np.sin(np.deg2rad(28)), -1.0, 703 * np.cos(np.deg2rad(28))],
    0,
)  # Abductor muscle force with 35 deg


# print('Total applied Hip-joint contact force = ', np.sqrt(np.square(-2317*np.sin(np.deg2rad(30))) + np.square(-2317*np.cos(np.deg2rad(30)))))
# print('Total applied Abductor muscle force = ', np.sqrt(np.square(-703*np.sin(np.deg2rad(35))) + np.square(703*np.cos(np.deg2rad(35)))))

nEl = data.nael
# rigidVol = data.nrel * 10.0
solidVol = data.nsel * 1.0
materialvolumefraction = 0.5

# Calculate the objective function
# objective input: (design variable value, SED)


def objective(comp, sumXp, volfrac):
    return comp


def sensitivity(xp, uKu, penal):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu


def constraint(comp, sumXp, volfrac):
    return (sumXp - solidVol) / nEl - materialvolumefraction
    # return 0.0


def constraintSensitivity(xp, uKu, penal):
    return 1.0 / nEl
    # return 0.0


# Callback implementation
data.obj(objective)
data.objsens(sensitivity)

# Define constraint
data.cons(constraint)
data.conssens(constraintSensitivity)


# Use local volume constraint additionally
# Local volume constraint input: (Rlocvol, alpha)
# data.localVolume(0.5, 0.6)

# Homogeneous initial condition
data.initialcondition(materialvolumefraction)

# Output vtr
data.vtr(20)

# step 3:
# solve topopt problem with input data and wait for "complete" signal
complete = data.solve()
