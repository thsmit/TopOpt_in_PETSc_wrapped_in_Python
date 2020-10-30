#!/bin/bash

# run with problem name as input var
# for example run:
# ./run_topopt.sh sphere_full.py

### SETUP ENVIRONMENT
# by hand: env2lmod
module load gcc/4.8.5 cmake/3.16.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4

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
pwd

cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib.so .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/examples/$1 .

# ADJUSTABLE PARAMETERS
#EULER_MEMORY="4000"
#NCPU=32
#WALL_TIME="12:00"

EULER_MEMORY="2000"
NCPU=16
WALL_TIME="4:00"

# FUNCTION CALL
bsub -n ${NCPU} -W ${WALL_TIME} -R ib -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python $1

# VIEW JOBS
bjobs