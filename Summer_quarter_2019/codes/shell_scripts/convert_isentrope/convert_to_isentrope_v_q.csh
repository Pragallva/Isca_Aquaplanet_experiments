#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

for i in {0..5};
do       
python $sum19/codes/python_scripts/isca/convert_to_isentropic_coord/convert_to_isentropic_coord.py $1 0 'vcomp' 'sphum'
done;        

