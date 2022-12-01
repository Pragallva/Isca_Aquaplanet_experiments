#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

for i in {0..9};
do       
python $win19/codes/python_scripts/Extract_ncfile_save_fluxes_radiation_depths.py $1 $2 $i
done;        

python $win19/codes/python_scripts/Reload_save_interpolated_depths.py $1 $2 10
python $win19/codes/python_scripts/Reload_save_thermo_dynamic_mmc.py $1 $2 10

##python $sum18/codes/python_scripts/Reload_sfc_temp_save_interpolated_depths.py  aqua isca$1 10


 
### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
