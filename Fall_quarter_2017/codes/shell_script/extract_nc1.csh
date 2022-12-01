#!/bin/csh -f

#SBATCH --time=05:00:00
#SBATCH --partition=bigmem
#SBATCH --constraint=256G
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-user=pragallva@uchicago.edu
#SBATCH --mail-type=ALL

python $project2/codes/python_scripts/Extract_ncfile_save_fluxes_radiation.py $1 $2
python $project2/codes/python_scripts/Reload_save_interpolated.py $1 $2
python $project2/codes/python_scripts/plot_and_save_div_MSE.py $1 $2
python $project2/codes/python_scripts/interpolation_from_sigma_to_pres_coord.py $1 $2
python $project2/codes/python_scripts/plot_stream_wind_etc.py $1 $2

### This piece of code loads nc file, saves MSE fluxes, reloads the flux data, interpolates it and plots it for zonal mean. It also saves changes to pressure coordinates and saves  stream function. Also plots zonal mean and stream function. 
