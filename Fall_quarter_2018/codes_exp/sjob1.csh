#!/bin/csh -f

#SBATCH --time=35:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=6000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python aqua_const_depth.py $1

