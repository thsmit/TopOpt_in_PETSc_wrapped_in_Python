#!/usr/bin/env python3

# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

import topoptlib

data = topoptlib.Data()
data.structuredGrid(
    (0.0, 192.0, 0.0, 64.0, 0.0, 104.0, 1.0, 7.0, 0.0, 0.0, 0.0, 0.0), (193, 65, 105)
)
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

if data.nael == 414813 and data.nsel == 4175 and data.nrel == 10521:
    open("test_bracket SUCCESFULL", "w+")
else:
    open("test_bracket not succesfull!", "w+")
