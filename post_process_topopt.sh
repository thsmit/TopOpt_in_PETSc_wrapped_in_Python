#!/bin/bash

# Author: Thijs Smit, Dec 2020
# Copyright (C) 2020 ETH Zurich

# Disclaimer:
# The authors reserves all rights but does not guaranty that the code is
# free from errors. Furthermore, we shall not be liable in any event
# caused by the use of the program.

# Import module
module purge
module load python/2.7.6
module list

cd $1
pwd
mkdir vtu

python $HOME/TopOpt_in_PETSc/bin2vtu.py 0

# User input
echo Final dataset number? ...
read fin

# User input
echo Convert all or only final 1/0?
read input

# compile
if [ $input -eq 1 ]
    then
        for i in 00 01 02 03 04 05 06 07 08 09 10 11 20 30 40 50 60 70 80 90 100
        do
        python $HOME/TopOpt_in_PETSc/bin2vtu.py $i
        mv output_000$i.vtu vtu
        done
fi

python $HOME/TopOpt_in_PETSc/bin2vtu.py $fin
mv output_000$fin.vtu vtu
