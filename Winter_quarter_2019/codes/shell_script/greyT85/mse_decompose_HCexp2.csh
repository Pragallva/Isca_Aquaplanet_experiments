#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --nodes=4
#SBATCH --mem=128000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

for i in {0..9};
do
python $win19/codes/python_scripts/greyT85/Extract_ncfile_save_fluxes_radiation_depths.py $1 $i 720 1440
done;

##python $win19/codes/python_scripts/grey/Reload_save_interpolated_depths.py $1 10
##python $win19/codes/python_scripts/grey/Reload_save_thermo_dynamic_mmc.py $1 10

##python $sum18/codes/python_scripts/grey/Reload_sfc_temp_save_interpolated_depths.py  aqua isca$1 10


 
### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
