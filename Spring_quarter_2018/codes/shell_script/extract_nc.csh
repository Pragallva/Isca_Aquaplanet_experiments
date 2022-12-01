#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

#for i in {29};
#do 
#python $spring18/codes/python_scripts/Extract_ncfile_save_fluxes_radiation.py $1 $2 29  
# done;
python $spring18/codes/python_scripts/Reload_save_interpolated.py $1 $2 30
python $spring18/codes/python_scripts/plot_and_save_div_MSE.py $1 $2
#python $spring18/codes/python_scripts/interpolation_from_sigma_to_pres_coord.py $1 $2
#python $spring18/codes/python_scripts/plot_stream_wind_etc.py $1 $2

### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
