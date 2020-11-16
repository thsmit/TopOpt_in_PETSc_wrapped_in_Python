#!/usr/bin/env python3

# Author: Thijs Smit, June 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

import topoptlib
import numpy as np

data = topoptlib.Data()
data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0), (129, 65, 65))
#data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (161, 81, 81))
#data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0), (257, 129, 129))
Emin, Emax, nu, dens, penal = 1.0e-9, 1.0, 0.3, 1.0, 3.0
data.material(Emin, Emax, nu, dens, penal)
data.filter(2, 0.08)
data.mma(40)
data.loadcases(1)
data.bc(0, 1, [0, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
data.bc(0, 2, [0, 1, 2, 4], [2], [-0.001], 0)
data.bc(0, 2, [0, 1, 1, 2, 2, 4], [2], [-0.0005], 0)
data.bc(0, 2, [0, 1, 1, 3, 2, 4], [2], [-0.0005], 0)

materialvolumefraction = 0.12
nEl = data.nElements

def objective(comp, sumXp, xp, uKu):
    return comp

def sensitivity(sumXp, xp, uKu):
    return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu

def constraint(comp, sumXp, xp, uKu):
    return sumXp / nEl - materialvolumefraction

def constraintSensitivity(sumXp, xp, uKu):
    return 1.0 / nEl

data.obj(objective)
data.objsens(sensitivity)
data.cons(constraint)
data.conssens(constraintSensitivity)
data.initialcondition(materialvolumefraction)
complete = data.solve()

# get a function to run for testing
def Test(trueFX, runtime):
    trueFX = int(trueFX*1000)/1000
    runtime = int(runtime*1000)/1000
    if trueFX == 2.164:
        print(trueFX, runtime)
        open("test_beam SUCCESFULL","w+")
    else:
        open("test_beam not succesfull!","w+")

data.check(Test)

print('True Compliance should be: ', 2.164643, 3.361824, 8.570609)