#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --nodes=4
#SBATCH --mem=128000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

for i in {0..9};
do       
python $win19/codes/python_scripts/greyT85/lower_filter/Calculate_EKE_depths.py $1 $i 0 720
done;        
###python $win19/codes/python_scripts/grey/lower_filter/Reload_EKE_save_interpolated_depths.py $1 10

 
### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
