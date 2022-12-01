#!/bin/csh -f

#SBATCH --time=35:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

##python fixed_sst.py $1
##python aqua_depths.py $1
python aqua_const_depth.py  $1  #qf_amp depth
