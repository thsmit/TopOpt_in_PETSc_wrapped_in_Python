#!/bin/bash

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

for i in 00 01 02 03 04 05 06 07 08 09 10 11 20 30 40 50 60 70 80 90 100 
#for i in {0..70}
do
python $HOME/TopOpt_in_PETSc/bin2vtu.py $i
mv output_000$i.vtu vtu
done

python $HOME/TopOpt_in_PETSc/bin2vtu.py $fin
mv output_000$fin.vtu vtu
