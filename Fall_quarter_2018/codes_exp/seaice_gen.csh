#!/bin/csh -f

#SBATCH --time=35:00:00
#SBATCH --account=rossby
#SBATCH --partition=rossby
#SBATCH --constraint=ib
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python symmetric_SH.py $1

