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
    (0.0, 1.0, 0.0, 1.2, 0.0, 1.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), (161, 193, 193)
)
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


if data.nael == 5816564 and data.nvel == 74978 and data.nrel == 6698:
    open("test_sphere SUCCESFULL", "w+")
else:
    open("test_sphere not succesfull!", "w+")
