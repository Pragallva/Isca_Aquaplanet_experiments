#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

##for i in {55..59};
##do  

set -e
     
python $fall19/codes/python_scripts/echam/convert_to_isentropic_coord/convert_to_isentropic_coord1.py $1 $2 'u'
python $fall19/codes/python_scripts/echam/convert_to_isentropic_coord/convert_to_isentropic_coord1.py $1 $2 'v'
python $fall19/codes/python_scripts/echam/convert_to_isentropic_coord/convert_to_isentropic_coord1.py $1 $2 'omega'
python $fall19/codes/python_scripts/echam/convert_to_isentropic_coord/convert_to_isentropic_coord1.py $1 $2 'q'
python $fall19/codes/python_scripts/echam/convert_to_isentropic_coord/convert_to_isentropic_coord1.py $1 $2 'geopoth'
python $fall19/codes/python_scripts/echam/convert_to_isentropic_coord/convert_to_isentropic_coord1.py $1 $2 't'
##done;        
