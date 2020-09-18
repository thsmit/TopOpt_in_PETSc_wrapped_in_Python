.. summary-start

TopOpt_in_PETSc_wrapped_in_Python
===============

``TopOpt_in_PETSc_wrapped_in_Python`` provides a Python wrapper and extends the functionality of the TopOpt_in_PETSc framework [1]_.

.. summary-end

**WARNING: TopOpt_in_PETSc_wrapped_in_Python is still in a development stage**

.. not-in-documentation-start

Implemented functionality
--------


Examples
--------

- Cantilever beam in ``beam.py``
- Multi-loads in ``multiloads.py``   
- Roof support in ``roof.py``
- Torsion ball in ``sphere.py``
- The Jet bracket in ``bracket.py``
- Local volume constraint on the bracket example in ``local.py``

Installation
------------


Tests
------------

Implemented tests in ``/tests``:

- Testing standard MBB problem with maxItr of 40 ``test_beam.py``
- Testing the standard MBB problem with two line loads ``test_multiload.py``
- Testing continuation of penalization ``test_continuation.py``
- Testing heavyside projection filtering ``test_projection.py``
- Testing stl readin of design domain, rigid domain ``test_sphere.py``

On ETH Euler: use the `test_topopt.sh` for automated building and running the tests

Learner
--------



Running on ETH Euler
--------

.. code:: bash

    env2lmod
    module load gcc/4.8.5 cmake/3.16.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4
    cd TopOpt_in_PETSc_wrapped_in_Python
    mkdir build
    cd build
    cmake ..
    make
    cd ..
    bsub -n 8 mpirun -n 8 python bracket.py

Or use the ``run_topopt.sh`` for automated building and running
    

Citing 
--------

.. code:: bash

    ...


Origional code
--------

.. [1]

    Aage, N., Andreassen, E., & Lazarov, B. S. (2015). Topology optimization using PETSc: An easy-to-use, fully parallel, open source topology optimization framework. 
    Structural and Multidisciplinary Optimization, 51(3), 565–572. https://doi.org/10.1007/s00158-014-1157-0
