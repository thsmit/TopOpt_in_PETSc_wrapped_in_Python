.. summary-start

TopOpt_in_PETSc_wrapped_in_Python
===============

``TopOpt_in_PETSc_wrapped_in_Python`` provides a Python wrapper called ``topoptlib`` and extends the functionality of the TopOpt_in_PETSc framework [1]_ [2]_ [3]_.

The Python interface simplifies the problem definition, is expanding the potential user-base and facilitates the use of a large-scale topology optimization framework for educational purposes. Furthermore, the functionality of the topology optimization framework is extended which contributes to its usability to real-world design applications. The functionality is demonstrated via the cantilever beam, bracket- and torsion ball examples. Several tests are provided which can be used to verify the proper working and compare the performance of the user’s system setup.

.. summary-end

|pic1| |pic2|

.. |pic1| image:: img/bracket_crop.gif
    :width: 45%

.. |pic2| image:: img/michell_crop.gif
    :width: 45%

.. not-in-documentation-start

Implemented functionality
--------

Large scale, high-resolution topology optimization including:

- STL (file format for storing surface geometry) input files to define the design domain, solid-, void- and rigid regions and voxelization
- Exclusion of passive elements from the optimization
- Application of loads and constraints using parametrization functions
- Multi-load cases and multi-constraints
- User defined objective- and constraint functions
- Local-volume constraint for 3D printing infill and bone like structures
- Continuation strategy for the penalization value
- Robust formulation using three-field density projection
- Test scripts for code verification

Installation
------------

The framework should be compiled once, on a cluster or a desktop computer. A problem file can use the functionality of the framework without compiling thereafter. A Linux system is recommended. A Windows machine should also work, however not tested.
The framework uses [CMake](https://cmake.org) to compile. The following third party libraries are required and located using CMake's ``find_package``.

- [Git] (https://git-scm.com/)
- [PETSc](https://www.mcs.anl.gov/petsc/mirror/release-snapshots/): version 3.14.1
- [Python] (https://www.python.org/): version 3.7

It needs PETSc to be installed:

.. code:: bash

    # Download PETSc source from its `release-snapshots <https://www.mcs.anl.gov/petsc/mirror/release-snapshots/>`_.
    # Then follow the guide on `Quickest Quick-start <https://petsc.org/release/install/install_tutorial/#qqtw-quickest-quick-start-in-the-west>`_.
    cd petsc-3.14.1
    ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=0 --download-f2cblaslapack=1 --with-debugging=0
    # PETSC_DIR is the path to the PETSc just downloaded
    make PETSC_DIR=/home/ts/Documents/petsc PETSC_ARCH=arch-linux-c-opt all
    make PETSC_DIR=/home/ts/Documents/petsc PETSC_ARCH=arch-linux-c-opt check

To download this framework:

.. code:: bash

    git clone https://github.com/thsmit/TopOpt_in_PETSc_wrapped_in_Python.git

To compile the framework (paths will differ):

.. code:: bash

    export PETSC_ARCH=arch-linux-c-opt
    export PETSC_DIR=/home/ts/Documents/petsc
    cd TopOpt_in_PETSc_wrapped_in_Python
    mkdir build && cd build
    cmake .. -D PETSC_EXECUTABLE_RUNS=ON
    make -j $(nproc)

Running 'hello world' example
--------

Running a 'hello world' example from the command line. Generates standard cantilever beam and output .vtr files for viewing in Paraview.

.. code:: bash

    import topoptlib
    data = topoptlib.Data()
    data.solve()


Running examples
--------

To run the cantilever beam example on one CPU (adjust the problem's mesh according to the number of available CPU's):

.. code:: bash

    cd TopOpt_in_PETSc_wrapped_in_Python
    cp examples/beam.py .
    python3 beam.py

Available examples:

- Cantilever beam in ``beam.py``
- Multi-loads in ``multiloads.py``
- Torsion ball in ``sphere.py``
- The Jet engine bracket in ``bracket.py``


Running on ETH Euler (without installing PETSc)
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
    bsub -n 8 mpirun -n 8 python3 bracket.py

Or use ``run_topopt.sh`` for automated building and running


Tests
------------

Several tests are provided to verify the proper working of the framework. To run a test using 4 CPU's use:

.. code:: bash

    cd TopOpt_in_PETSc_wrapped_in_Python
    cp topoptlib/test/test_beam.py .
    mpirun -n 4 python3 test_beam.py

Implemented tests in ``/tests``:

- Testing standard MBB problem with maxItr of 40 ``test_beam.py``
- Testing the standard MBB problem with two line loads ``test_multiload.py``
- Testing continuation of penalization ``test_continuation.py``
- Testing heavyside projection filtering ``test_projection.py``
- Testing stl readin of design domain, rigid domain ``test_sphere.py``
- Testing stl readin of design domain, rigid domain ``test_bracket.py``
- Testing the robust approach ``test_michell.py``

Or use ``test_topopt.sh`` for automated building and running the tests

Post-processing (easy)
--------

The framework can write .vtr files of the designs with in point data. The designs can be viewed in Paraview (https://www.paraview.org/). The point data can be transformed into cell data by using Paraview's PointToCellData filter.
To generate .vtr files add the following command to the problem definition:

.. code:: bash

    vtr(20)


Post-processing (original)
--------

After solving the problem the output is written to a ``output.dat`` file. The designs can be viewed in Paraview (https://www.paraview.org/).
To generate .vtu files from the output file use ``post_process_topopt.sh`` with Python 2 (with * the file path and name where the output file is stored):

.. code:: bash

    cd TopOpt_in_PETSc_wrapped_in_Python
    ./post_process_topopt.sh *

Disclaimer
--------

The authors reserves all rights but does not guaranty that the code is free from errors. Furthermore, we shall not be liable in any event caused by the use of the program.

Citing
--------

For citing this work use:

.. code:: bib

    @article{Smit2021,
    author = {Smit, Thijs and Aage, Niels and Ferguson, Stephen J and Helgason, Benedikt},
    title = {{Topology optimization using PETSc : a Python wrapper and extended functionality}},
    journal = {Structural and Multidisciplinary Optimization},
    year = {2021}
    publisher = {Springer Berlin Heidelberg},
    doi = {10.1007/s00158-021-03018-7},
    url = {https://doi.org/10.1007/s00158-021-03018-7},
    }


Original code
--------

.. [1]

    Aage, N., Andreassen, E., & Lazarov, B. S. (2015). Topology optimization using PETSc: An easy-to-use, fully parallel, open source topology optimization framework.
    Structural and Multidisciplinary Optimization, 51(3), 565–572. https://doi.org/10.1007/s00158-014-1157-0

.. [2]

    http://topopt.dtu.dk/PETSc

.. [3]

    https://github.com/topopt/TopOpt_in_PETSc
