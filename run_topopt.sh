#!/bin/bash

# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

# run with problem name as input var
# for example run:
# ./run_topopt.sh sphere_full.py

### SETUP ENVIRONMENT
# by hand: env2lmod
#module load gcc/4.8.5 cmake/3.16.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4
module load gcc/4.8.5 cmake/3.16.5 openmpi/3.0.1 python/3.7.4

export PETSC_ARCH=arch-linux-c-opt
export PETSC_DIR=$home/petsc

# Print loaded modules
module list

# User input
echo Press 1 for compiling...
read var

# compile
if [ $var -eq 1 ]
    then 
        rm -rf topoptlib.so
        rm -rf build
        mkdir build
        cd build
        cmake ..
        make
        cd ..
fi

# RUN
cd $SCRATCH/wrapped;

id=`date '+%Y%m%d_%H:%M:%S'`;
#echo $id

mkdir $id
cd $id
#pwd

cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib.so .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/examples/$1 .

# ADJUSTABLE PARAMETERS
EULER_MEMORY="4000"
NCPU=32
WALL_TIME="12:00"

#EULER_MEMORY="1000"
#NCPU=4
#WALL_TIME="01:00"

# FUNCTION CALL
bsub -n ${NCPU} -W ${WALL_TIME} -R ib -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python $1

# VIEW JOBS
bjobs
