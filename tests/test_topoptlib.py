#!/usr/bin/python3

# Author: Thijs Smit, May 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

def test_topoptlib():

    import topoptlib

    # run default problem
    data = topoptlib.Data()
    data.structuredGrid((0.0, 2.0, 0.0, 1.0, 0.0, 1.0), (65, 33, 33))
    Emin = 1.0e-9
    Emax = 1.0
    penal = 3.0
    data.material(Emin, Emax, 0.3, penal)
    data.filter(1, 0.08)
    data.mma(4)
    complete = data.solve()

    assert complete == 1, "Problem not completed"
