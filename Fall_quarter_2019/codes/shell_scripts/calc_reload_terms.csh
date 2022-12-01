#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --partition=tas1
#SBATCH --constraint=ib
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python $sum19/codes/python_scripts/isca/Reload_save_interpolated_vdecomp_momentum.py 10x_stress 10 30
python $sum19/codes/python_scripts/isca/Reload_save_interpolated_depths.py 10x_stress 10 30

python $sum19/codes/python_scripts/isca/Reload_save_interpolated_vdecomp_momentum.py 0.1x_stress 10 30
python $sum19/codes/python_scripts/isca/Reload_save_interpolated_depths.py 0.1x_stress 10 30

python $sum19/codes/python_scripts/isca/Reload_save_interpolated_vdecomp_momentum.py 2x_stress 10 30
python $sum19/codes/python_scripts/isca/Reload_save_interpolated_depths.py 2x_stress 10 30

python $sum19/codes/python_scripts/isca/Reload_save_interpolated_vdecomp_momentum.py 0.5x_stress 10 30
python $sum19/codes/python_scripts/isca/Reload_save_interpolated_depths.py 0.5x_stress 10 30

python $sum19/codes/python_scripts/isca/Reload_save_interpolated_vdecomp_momentum.py HC0_la5m_oc5m 10 30
python $sum19/codes/python_scripts/isca/Reload_save_interpolated_depths.py HC0_la5m_oc5m 10 30
