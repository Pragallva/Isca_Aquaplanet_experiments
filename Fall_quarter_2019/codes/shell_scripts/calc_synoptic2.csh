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
python $win19/codes/python_scripts/isca/calculate_synoptic_variability.py $1 $i 'high_pass' 30
done;        
python $win19/codes/python_scripts/isca/Reload_save_synoptic_variability.py $1 10 'high_pass' 30




