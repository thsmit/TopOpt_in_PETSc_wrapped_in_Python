#!/bin/bash

# Import module
module purge
module load python/2.7.6
module list

cd $1
pwd 
mkdir vtu

for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
do
python $HOME/TopOpt_in_PETSc/bin2vtu.py $i
mv output_000$i.vtu vtu
done
