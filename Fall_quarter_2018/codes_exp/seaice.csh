#!/bin/csh -f

#SBATCH --time=35:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python seaice_experiments.py land_ice 1 $1

