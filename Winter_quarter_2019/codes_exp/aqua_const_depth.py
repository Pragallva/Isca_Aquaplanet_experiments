import numpy as np
import os

import sys
dirc=sys.argv

from gfdl import Experiment, DiagTable

experiment_name=str(dirc[1])+'m'
baseexp = Experiment(experiment_name, overwrite_data=False)

baseexp.inputfiles = ['/project/tas1/pragallva/Fall_quarter_2017/Isca/input/ozone_1990.nc'] 

diag_spinup = DiagTable()
diag_spinup.add_file('atmos_4xdaily', 6, 'hours', time_units='hours')

diag_spinup.add_field('dynamics', 'ps', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'pres_full', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'pres_half', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'bk', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'pk', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'ucomp', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'vcomp', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'temp', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'sphum', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'omega', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'zsurf', time_avg=True, files=['atmos_4xdaily']) ## geopotential height at the surface
diag_spinup.add_field('dynamics', 'height', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'height_half', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'EKE', time_avg=True, files=['atmos_4xdaily'])  
diag_spinup.add_field('dynamics', 'vor', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('dynamics', 'div', time_avg=True, files=['atmos_4xdaily'])

diag_spinup.add_field('atmosphere', 'precipitation', time_avg=True, files=['atmos_4xdaily'])

diag_spinup.add_field('rrtm_radiation', 'flux_sw', time_avg=True, files=['atmos_4xdaily']) # Net SW surface flux, in W/m^2
diag_spinup.add_field('rrtm_radiation', 'flux_lw', time_avg=True, files=['atmos_4xdaily']) # Net LW surface flux, in W/m^2
diag_spinup.add_field('rrtm_radiation', 'rrtm_albedo', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('rrtm_radiation', 'olr', time_avg=True, files=['atmos_4xdaily'])   # OLR W/m^2,             in W/m^2
diag_spinup.add_field('rrtm_radiation', 'coszen', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('rrtm_radiation', 'toa_sw', time_avg=True, files=['atmos_4xdaily']) # Net TOA SW flux, W/m2 
diag_spinup.add_field('mixed_layer', 'flux_oceanq', time_avg=True, files=['atmos_4xdaily']) ## Ocean heat flux in Watts/m^2
diag_spinup.add_field('mixed_layer', 't_surf', time_avg=True, files=['atmos_4xdaily'])  ## Surface temperature in K
diag_spinup.add_field('mixed_layer', 'flux_lhe', time_avg=True, files=['atmos_4xdaily'])  ## 'latent heat flux up at surface','watts/m2'
diag_spinup.add_field('mixed_layer', 'flux_t', time_avg=True, files=['atmos_4xdaily']) ## sensible heat flux up at surface, 'W/m^2'
diag_spinup.add_field('mixed_layer', 'ml_heat_cap', time_avg=True, files=['atmos_4xdaily']) ## sensible heat flux up at surface, 'W/m^2'


baseexp.use_diag_table(diag_spinup)

baseexp.compile()

baseexp.clear_rundir()


baseexp.namelist['main_nml'] = {
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':360,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
}

baseexp.namelist['surface_flux_nml']['use_virtual_temp'] = True
baseexp.namelist['surface_flux_nml']['entr_ratio'] = 0.0
baseexp.namelist['surface_flux_nml']['parcel_buoy'] = 0.0
baseexp.namelist['spectral_dynamics_nml']['damping_coeff'] = 6.9444444e-05
baseexp.namelist['spectral_dynamics_nml']['robert_coeff']  = 0.04
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.1
baseexp.namelist['spectral_dynamics_nml']['scale_heights'] = 5.0
baseexp.namelist['spectral_dynamics_nml']['exponent'] = 2.0
baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 50


baseexp.namelist['mixed_layer_nml']['depth'] = float(dirc[1])
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25 #0.3
baseexp.namelist['spectral_dynamics_nml']['ocean_topog_smoothing'] = 0.8
baseexp.namelist['spectral_dynamics_nml']['raw_filter_coeff'] = 1.
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150.
baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 3600
baseexp.namelist['rrtm_radiation_nml']['do_read_ozone']=True
baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False

baseexp.namelist['mixed_layer_nml']['load_qflux'] = False
baseexp.namelist['mixed_layer_nml']['do_qflux'] = False

baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True
baseexp.namelist['idealized_moist_phys_nml']['turb'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_damping'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_virtual'] = False
baseexp.namelist['idealized_moist_phys_nml']['roughness_mom'] = 2e-04
baseexp.namelist['idealized_moist_phys_nml']['roughness_heat'] = 2e-04
baseexp.namelist['idealized_moist_phys_nml']['roughness_moist'] = 2e-04
baseexp.namelist['spectral_dynamics_nml']['lon_max'] = 256
baseexp.namelist['spectral_dynamics_nml']['lat_max'] = 128
baseexp.namelist['spectral_dynamics_nml']['num_fourier'] = 85
baseexp.namelist['spectral_dynamics_nml']['num_spherical'] = 86

##baseexp.namelist['qe_moist_convection_nml']['Tmin']=50
##baseexp.namelist['qe_moist_convection_nml']['Tmax']=450


## Radiative gases mixixng ratio
baseexp.namelist['rrtm_radiation_nml']['include_secondary_gases'] = True
baseexp.namelist['rrtm_radiation_nml']['temp_lower_limit'] = 50
baseexp.namelist['rrtm_radiation_nml']['temp_upper_limit'] = 450
baseexp.namelist['rrtm_radiation_nml']['co2ppmv'] = 348.0
baseexp.namelist['rrtm_radiation_nml']['ch4_val'] = 1650.e-9
baseexp.namelist['rrtm_radiation_nml']['n2o_val'] = 306e-9
baseexp.namelist['rrtm_radiation_nml']['cfc11_val'] = 124.7e-12
baseexp.namelist['rrtm_radiation_nml']['cfc12_val'] = 232.9e-12
baseexp.namelist['constants_nml']['omega'] = 7.2921150e-5

baseexp.screen_runmonth_prefix = 'aqua_'+str(dirc[1])+'m'

baseexp.runmonth(1, use_restart=False, num_cores=16)
for i in range(2, 361):  
    baseexp.runmonth(i, num_cores=16)
