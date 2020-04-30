.. summary-start

TopOpt_in_PETSc_wrapped_in_Python
===============

|Topology Optimization| |Python| |Large scale| |HPC| |Python wrapper| |PETSc| |Parralel computing|

``topoptlib`` provides a Python wrapper and extens the functionality of the TopOpt in PETSc framework.

.. summary-end

.. not-in-documentation-start

Examples
--------

- [ ] Cantilever beam in ``beam.py``
- [ ] Roof support in ``roof.py``
- [ ] Force inverter in ``inverter.py``


Implemented functionality
----------------------

- ``loadcases(.., .., ..)`` multipleloads
- ``passive()`` passive elements
- ``vtu()`` generate vtu files for paraview
- ``stl()`` generate stl file of final design for 3D-printing


Installation of topoptlib
------------

Requirements:
- PETSc 3.13.0
- Python 3.7.4

Tests
------------

Implemented tests in ``/tests``:
- [x] Testing standard MBB problem with maxItr of 4 ``test_topoptlib.py``
- [ ] Testing multiloads functionality
- [ ] Testing passive elements functionality

Running tests ``make test``


ToDo
--------

- [x] Multiloads
- [ ] Passive elements
- [ ] Multiple constraints
- [ ] Wrap BC
- [ ] Test with Petsc 3.13.0
- [ ] Update tests
- [ ] Python multithreading
- [ ] Default variable values
- [ ] vtu gen
- [ ] stl gen
- [ ] performance script
- [ ] infill constraint
- [ ] learner with jupyter notebook


Running on ETH Euler
--------

.. code:: bash

    env2lmod
    module load gcc/4.8.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4
    cd TopOpt_in_PETSc_wrapped_in_Python
    make topoptlib
    make test
    make myclean
    bsub -n 8 mpirun python problem_definition_file.py

Origional code
--------

.. code:: bash

    `Aage, N., Andreassen, E., & Lazarov, B. S. (2015). Topology optimization using PETSc: An easy-to-use, 
    fully parallel, open source topology optimization framework. Structural and Multidisciplinary Optimization, 51(3), 565–572. 
    https://doi.org/10.1007/s00158-014-1157-0

Citing 
--------

`...`
