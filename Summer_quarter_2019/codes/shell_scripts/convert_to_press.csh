#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

##for i in {0..9};
##do       
python $sum19/codes/python_scripts/am2/convert_to_pressure_coord/convert_to_pressure_coord.py $1 $2
##done;        

