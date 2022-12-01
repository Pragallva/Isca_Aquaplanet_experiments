#!/bin/csh -f

#SBATCH --time=5:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python $project2/codes/python_scripts/interpolation_from_sigma_to_pres_coord.py $1 $2
