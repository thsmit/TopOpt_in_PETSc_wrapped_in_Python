#!/usr/bin/env python3

# Author: Thijs Smit, Dec 2020
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
data.structuredGrid((0.0, 192.0, 0.0, 64.0, 0.0, 104.0, 1.0, 7.0, 0.0, 0.0, 0.0), (193, 65, 105))

# readin STL file in binary format
# stl read: (encoding, backround, treshold, box around stl: (min corner)(max corner), full path to file)
# Passive elements: 1.0
# Active elements: -1.0
# Solid elements: 2.0
# Rigid elements: 3.0
# Do not overwrite: 0.0
data.stlread(-1.0, 1.0, 8, (-23.0, -1.0, -103.0), (169.0, 63.0, 1.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineDesignDomainFine.stl')
data.stlread(2.0, 0.0, 4, (-23.0, -1.0, -103.0), (169.0, 63.0, 1.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineSolidDomainFine.stl')
data.stlread(3.0, 0.0, 8, (-23.0, -1.0, -103.0), (169.0, 63.0, 1.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineRigidDomainFine.stl')

# readin STL file in binary format
# stl read: ((box around stl: (min corner)(max corner))full path to file)
#data.stlread_domain((-23.0, -1.0, -103.0), (169.0, 63.0, 1.0), '/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineDesignDomainFine.stl')

# stl read: load the solid domain into the same coordinate system as the design domain
#data.stlread_solid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineSolidDomainFine.stl')

# stl read: load the solid domain into the same coordinate system as the design domain
#data.stlread_rigid('/cluster/home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/stl/jetEngineRigidDomainFine.stl')

# Optional printing:
#print(data.nNodes)
#print(data.nElements)
#print(data.nDOF)
#print(data.nael)
#print(data.nsel)
#print(data.nrel)


if data.nael == 414813 and data.nsel == 4175 and data.nrel == 10521:
    open("test_bracket SUCCESFULL","w+")
else:
    open("test_bracket not succesfull!","w+")