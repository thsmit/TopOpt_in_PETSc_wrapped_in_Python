#!/usr/bin/python3

# Author: Thijs Smit, May 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

def main():
    import topoptlib
    import numpy as np

    # step 1:
    # Create data class to store input data
    data = topoptlib.Data()

    # step 2:
    # define input data
    # mesh: (domain: x, y, z, center)(mesh: number of nodes)
    #data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0), (65, 33, 33))
    data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0), (129, 65, 65))
    #data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0), (257, 129, 129))
    
    # Optional printing:
    #print(data.nNodes)
    #print(data.nElements)
    #print(data.nDOF)

    # material: (Emin, Emax, nu, penal)
    Emin, Emax, nu, Dens, penal = 1.0e-6, 1.0, 0.3, 1.0, 3.0
    data.material(Emin, Emax, nu, Dens, penal)

    # setup continuation of penalization: (Pinitial, Pfinal, stepsize)
    #data.continuation()

    # setup heavyside projection filter (betaFinal, stepsize, eta)
    #data.projection() 

    # filter: (type, radius)
    # filter types: sensitivity = 0, density = 1, 
    # using 0.08, 0.04 or 0.02
    data.filter(1, 0.04)

    # optimizer: (maxIter)
    data.mma(400)

    # loadcases: (# of loadcases)
    data.loadcases(1)

    # bc: (loadcase, type, [checker: lcoorp[i+?], xc[?]], [setter: dof index], [setter: values])
    data.bc(0, 1, [0, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
    data.bc(0, 2, [0, 1, 2, 4], [2], [-0.001], 0)
    data.bc(0, 2, [0, 1, 1, 2, 2, 4], [2], [-0.0005], 0)
    data.bc(0, 2, [0, 1, 1, 3, 2, 4], [2], [-0.0005], 0)

    #data.bc(1, 1, [0, 0], [0, 1, 2], [0.0, 0.0, 0.0], 0)
    #data.bc(1, 2, [0, 1, 1, 3], [1], [0.001], 0)
    #data.bc(1, 2, [0, 1, 1, 3, 2, 4], [1], [0.0005], 0)
    #data.bc(1, 2, [0, 1, 1, 3, 2, 5], [1], [0.0005], 0)

    materialvolumefraction = 0.12
    #complicancetarget = 1.5
    nEl = data.nElements

    # Calculate the objective function, senitivity, constraint and constraint sensitivity
    def objective(comp, sumXp):
        return comp # for minimizing compliance
        #return sumXp # for minimizing volume

    def sensitivity(xp, uKu):
        return -1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu # for minimizing compliance
        #return 1.0 # for minimizing volume

    def constraint(comp, sumXp):
        #print('sumXp: ', sumXp)
        #return 0.0
        return sumXp / nEl - materialvolumefraction # for minimizing compliance
        #return comp / complicancetarget - 1.0 # for minimizing volume

    def constraintSensitivity(xp, uKu):
        #return 0.0
        return 1.0 / nEl # for minimizing compliance
        #return (-1.0 * penal * np.power(xp, (penal - 1)) * (Emax - Emin) * uKu) / complicancetarget # for minimizing volume

    # Callback implementatio
    data.obj(objective)
    data.objsens(sensitivity)

    # Define constraint
    data.cons(constraint)
    data.conssens(constraintSensitivity)

    # Use local volume constraint additionally
    # Local volume constraint input: (Rlocvol, alpha)
    #data.localVolume(0.16, 0.12)

    # Volume constraint is standard, input (volume fraction)
    data.initialcondition(materialvolumefraction)

    # number of cores used
    nc = 32

    # get a function to run for testing
    def Test(trueFX, runtime, memory):
        trueFX = int(trueFX*1000)/1000
        CPUtime = int (np.round(runtime * nc, decimals=0)) # in seconds
        memory= int( (memory*0.000001) * nc)  # in mega bytes
        print('trueFX: ', trueFX)
        print('CPU time: ', CPUtime, ' sec')
        print('Memory use per core, extrapolated: ', memory, ' MB')

    data.check(Test)

    # step 3:
    # solve topopt problem with input data and wait for "complete" signal
    complete = data.solve()

if __name__ == "__main__":
    main()