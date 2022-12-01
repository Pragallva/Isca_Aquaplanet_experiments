#!/bin/csh -f

#SBATCH --time=00:08:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@@uchicago.edu
#SBATCH --mail-type=ALL

python $project1/Isca/postprocessing/plevel_interpolation/scripts/run_plevel.py $1
