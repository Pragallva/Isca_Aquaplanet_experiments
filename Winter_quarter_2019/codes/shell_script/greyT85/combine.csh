#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --nodes=4
#SBATCH --mem=128000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL


for i in {0..9};
do
python $win19/codes/python_scripts/greyT85/combine_6_months.py $1 $i $2
done;    

 
