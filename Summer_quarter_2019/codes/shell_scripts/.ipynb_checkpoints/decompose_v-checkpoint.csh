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
##python $sum19/codes/python_scripts/isca/calc_eddy_fluxes/momentum_equation_terms.py $1 $i 30
##done;        
python $sum19/codes/python_scripts/isca/calc_eddy_fluxes/v_decomposition.py $1 10 30

