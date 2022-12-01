#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --nodes=4
#SBATCH --mem=128000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL


##for i in {0..9};
##do
##python $win19/codes/python_scripts/greyT85/combine_6_months.py $1 $i 1
##python $win19/codes/python_scripts/greyT85/combine_6_months.py $1 $i 0
##done;
       
##python $win19/codes/python_scripts/greyT85/Reload_save_interpolated_depths.py $1 10
##python $win19/codes/python_scripts/greyT85/Reload_save_thermo_dynamic_mmc.py  $1 10       
##python $win19/codes/python_scripts/greyT85/lower_filter/Reload_EKE_save_interpolated_depths.py  $1 10
python $win19/codes/python_scripts/greyT85/Calculate_baroclinicity.py  $1 10

