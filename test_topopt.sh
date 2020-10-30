#!/bin/bash


### SETUP ENVIRONMENT
# by hand: env2lmod
module load gcc/4.8.5 cmake/3.16.5 openmpi/3.0.1 petsc/3.10.5 python/3.7.4


# Print loaded modules
module list

# compile
rm -rf topoptlib.so
rm -rf build
mkdir build
cd build
cmake ..
make
cd ..

# RUN
cd $SCRATCH/wrapped;

id=`date '+%Y%m%d_%H:%M:%S_test'`;
#echo $id

mkdir $id
cd $id
#pwd

cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib.so .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_beam.py .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_beam_local.py .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_bracket.py .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_multiload.py .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_sphere.py .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_continuation.py .
cp ../../../../home/thsmit/TopOpt_in_PETSc_wrapped_in_Python/topoptlib/test/test_projection.py .

# ADJUSTABLE PARAMETERS
EULER_MEMORY="1000"
NCPU=8
WALL_TIME="1:00"

# FUNCTION CALL
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_beam.py
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_beam_local.py
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_bracket.py
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_multiload.py
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_sphere.py
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_continuation.py
bsub -n ${NCPU} -W ${WALL_TIME} -R "rusage[mem=${EULER_MEMORY}]" mpirun -n ${NCPU} python test_projection.py

# VIEW JOBS
bjobs