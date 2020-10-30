#!/bin/bash

# Import module
module purge
module load python/2.7.6
module list

cd $1
pwd 
mkdir vtu

for i in {0..40}
do
python $HOME/TopOpt_in_PETSc/bin2vtu.py $i
mv output_000$i.vtu vtu
done
