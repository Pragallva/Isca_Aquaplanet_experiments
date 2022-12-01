#!/bin/csh -f

#SBATCH --time=35:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16 
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python HC_experiments.py $1 $2 $3

## land_depth, ocean_depth and edge
