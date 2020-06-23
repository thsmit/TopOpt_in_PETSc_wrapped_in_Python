.. summary-start

TopOpt_in_PETSc_wrapped_in_Python
===============

``TopOpt_in_PETSc_wrapped_in_Python`` provides a Python wrapper and extends the functionality of the TopOpt_in_PETSc framework.

.. summary-end

**WARNING: adaptive is still in a beta development stage**

.. not-in-documentation-start

Examples
--------

- Cantilever beam in ``beam.py``
- Multi-loads in ``multiloads.py``
<img src="![3](https://user-images.githubusercontent.com/52911749/85385983-48616480-b543-11ea-86a3-45c1efe5d447.png)"> </img>
- Roof support in ``roof.py``

Installation
------------


Tests
------------

Implemented tests in ``/tests``:

- Testing standard MBB problem with maxItr of 40 ``test_beam.py``

Running tests:

.. code:: bash

    module load gcc/4.8.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4
    cd TopOpt_in_PETSc_wrapped_in_Python
    make topoptlib
    make test
    make myclean


Learner
--------


Wish-list
--------

- [ ] STL input
- [ ] Passive elements proper implementation
- [ ] continuation strategy
- [ ] infill constraint
- [ ] Automatic STL and VTK output
- [ ] performance script
- [ ] Multiple constraints
- [ ] Test with Petsc 3.13.0

Running on ETH Euler
--------

.. code:: bash

    env2lmod
    module load gcc/4.8.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4
    cd TopOpt_in_PETSc_wrapped_in_Python
    make topoptlib
    make test
    bsub -n 8 mpirun -n 8 python multiloads.py
    make myclean

Origional code
--------

.. code:: bash

    Aage, N., Andreassen, E., & Lazarov, B. S. (2015). Topology optimization using PETSc: An easy-to-use, fully parallel, open source topology optimization framework. 
    Structural and Multidisciplinary Optimization, 51(3), 565–572. https://doi.org/10.1007/s00158-014-1157-0

Citing 
--------

.. code:: bash

    ...
