# Limitations: work on Linux machine
# Have available in your environment: GCC 4.8.5, OpenMPI 4.0.2, Python 3.7.4

#PETSC_ARCH=arch-linux-c-debug
#PETSC_DIR=/cluster/home/thsmit/petsc

NUMPY_INCLUDE = /cluster/apps/nss/python/3.7.4/x86_64/lib64/python3.7/site-packages/numpy/core/include/numpy

CFLAGS=-I. -I${PYTHON_ROOT}/include/python3.7m -I${NUMPY_INCLUDE}
FFLAGS=
CPPFLAGS=-I. -I${PYTHON_ROOT}/include/python3.7m -I${NUMPY_INCLUDE}
FPPFLAGS=
LOCDIR= 
EXAMPLESC=
EXAMPLESF=
MANSEC=
CLEANFILES=
NP=

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

topoptlib: wrapper.o loop.o TopOpt.o LinearElasticity.o MMA.o Filter.o PDEFilter.o MPIIO.o chkopts
	rm -rf topoptlib.so
	mpicc -shared -fPIC -o topoptlib.so wrapper.o loop.o TopOpt.o LinearElasticity.o MMA.o Filter.o PDEFilter.o MPIIO.o ${PETSC_SYS_LIB}
	${RM} wrapper.o loop.o TopOpt.o LinearElasticity.o MMA.o Filter.o PDEFilter.o MPIIO.o 
	rm -rf *.o

myclean:
	rm -rf main *.a *.so *.o output* binary* log* makevtu.pyc Restart*

test:
	bsub -I -n 8 mpirun pytest
	
