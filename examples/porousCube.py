#!/usr/bin/env python3

# Author: Yijun Zhou, Aug 2021
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.


def main():
    import numpy as np

    import topoptlib

    # step 1:
    # Create data class to store input data
    data = topoptlib.Data()

    # step 2:
    # define input data
    # mesh: (domain: x, y, z)(mesh: number of nodes)
    # data.structuredGrid((0.0, 96.0, 0.0, 60.0, 0.0, 96.0), (321, 201, 321))
    # data.structuredGrid((-25.0, 25.0, 25.0, 25.0, 0.0, 50.0), (65, 65, 65))
    data.structuredGrid((0.0, 5.0, 0.0, 5.0, 0.0, 5.0), (65, 65, 65))

    # readin STL file in binary format
    # stl read: (encoding, backround, threshold, box around stl: (min corner)(max corner), full path to file)
    # Passive elements: 1.0
    # Active elements: -1.0
    # Solid elements: 2.0
    # Rigid elements: 3.0
    # Void elements: 4.0
    # Initial condition: 5.0
    # Do not overwrite: 0.0
    data.stlread(
        5.0,
        -1.0,
        8,
        (-2.5, -2.5, 0.0),
        (2.5, 2.5, 5.0),
        "/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/porousCube.stl",
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
    # data.projection(64.0, 1.0, 0.5)
    # # robust formulation (betaFinal, betaInit, eta, delta)
    # data.robust(64.0, 1.0, 0.5, 0.1)

    # filter: (type, radius)
    # filter types: sensitivity = 0, density = 1, pde = 2,
    data.filter(1, 0.1)

    # optimizer: (maxIter, tol)
    data.mma(20, 0.01)

    # loadcases: (# of loadcases)
    data.loadcases(1)

    # bc: (loadcase, type, [coordinate axis], [coordinate value], [coordinate axis], [bc value], parametrization)
    # data.bc(0, 1, [2], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
    # fully constraint at mid-diaphysis
    # data.bc(0, 2, [0, 1, 2], [0, 0, 45.0], [2], [-1.0], 0)

    data.bc(0, 1, [2], [0.0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
    data.bc(0, 2, [2], [5.0], [2], [-0.001], 0)
    data.bc(0, 2, [2, 1], [5.0, 0.0], [2], [-0.0005], 0)
    data.bc(0, 2, [2, 1], [5.0, 5.0], [2], [-0.0005], 0)
    data.bc(0, 2, [2, 0], [5.0, 0.0], [2], [-0.0005], 0)
    data.bc(0, 2, [2, 0], [5.0, 5.0], [2], [-0.0005], 0)
    data.bc(0, 2, [2, 1, 0], [5.0, 0.0, 0.0], [2], [-0.00025], 0)
    data.bc(0, 2, [2, 1, 0], [5.0, 0.0, 5.0], [2], [-0.00025], 0)
    data.bc(0, 2, [2, 1, 0], [5.0, 5.0, 0.0], [2], [-0.00025], 0)
    data.bc(0, 2, [2, 1, 0], [5.0, 5.0, 5.0], [2], [-0.00025], 0)

    materialvolumefraction = 0.12
    nEl = data.nael
    # rigidVol = data.nrel * 10.0
    # solidVol = data.nsel * 1.0

    # Calculate the objective function
    # objective input: (design variable value, SED)

    def objective(comp, sumXp, volfrac):
        return comp

    def sensitivity(xp, uKu, penal):
        return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

    def constraint(comp, sumXp, volfrac):
        return sumXp / nEl - materialvolumefraction
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
    # data.localVolume(2.5, 0.6)

    # Homogeneous initial condition
    data.initialcondition(0.5)

    # Output vtr
    data.vtr(20)

    # step 3:
    # solve topopt problem with input data and wait for "complete" signal
    # complete = data.solve()
    data.solve()


if __name__ == "__main__":
    main()
