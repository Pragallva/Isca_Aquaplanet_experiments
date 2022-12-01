#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --nodes=4
#SBATCH --mem=128000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

for i in {0..9};
do       
python $win19/codes/python_scripts/greyT85/convert_to_pressure_coord.py $1 $i 720 1440
done;        
##python $win19/codes/python_scripts/grey/Calculate_baroclinicity.py $1 10

 
### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
