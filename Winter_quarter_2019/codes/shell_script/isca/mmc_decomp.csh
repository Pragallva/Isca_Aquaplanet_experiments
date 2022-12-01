#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

##for i in {0..9};
##do       
##python $win19/codes/python_scripts/isca/thermo_dynamic_decompose.py $1 $i
##done;        
python $win19/codes/python_scripts/isca/Reload_save_thermo_dynamic_mmc.py $1 10

 
### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
