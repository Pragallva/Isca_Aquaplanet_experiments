import numpy as np
import os

import sys
dirc=sys.argv

from gfdl import Experiment, DiagTable


qf=dirc[1]
depth=dirc[2]

experiment_name = "Lenka_oc"+str(depth)+"qf"+str(qf)+"T85"
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
diag_spinup.add_field('two_stream', 'swdn_sfc', time_avg=True, files=['atmos_4xdaily']) # Absorbed SW at surface, in W/m^2
diag_spinup.add_field('two_stream', 'swdn_toa', time_avg=True, files=['atmos_4xdaily']) # SW flux down at TOA, in W/m^2
diag_spinup.add_field('two_stream', 'lwup_sfc', time_avg=True, files=['atmos_4xdaily']) # LW flux up at surface, in W/m^2
diag_spinup.add_field('two_stream', 'lwdn_sfc', time_avg=True, files=['atmos_4xdaily']) # LW flux down at surface, in W/m^2
diag_spinup.add_field('two_stream', 'net_lw_sfc', time_avg=True, files=['atmos_4xdaily']) # Net upward LW flux at surface, in W/m^2
diag_spinup.add_field('two_stream', 'flux_lw', time_avg=True, files=['atmos_4xdaily']) # Net longwave radiative flux (positive up), in W/m^2
diag_spinup.add_field('two_stream', 'flux_sw', time_avg=True, files=['atmos_4xdaily']) # Net shortwave radiative flux (positive up), in W/m^2
diag_spinup.add_field('two_stream', 'olr', time_avg=True, files=['atmos_4xdaily']) # Net LW surface flux, in W/m^2
diag_spinup.add_field('two_stream', 'flux_rad', time_avg=True, files=['atmos_4xdaily']) # Total radiative flux (positive up), in W/m^2
diag_spinup.add_field('mixed_layer', 'flux_oceanq', time_avg=True, files=['atmos_4xdaily']) ## Ocean heat flux in Watts/m^2
diag_spinup.add_field('mixed_layer', 't_surf', time_avg=True, files=['atmos_4xdaily'])  ## Surface temperature in K
diag_spinup.add_field('mixed_layer', 'flux_lhe', time_avg=True, files=['atmos_4xdaily'])  ## 'latent heat flux up at surface','watts/m2'
diag_spinup.add_field('mixed_layer', 'flux_t', time_avg=True, files=['atmos_4xdaily']) ## sensible heat flux up at surface, 'W/m^2'
diag_spinup.add_field('mixed_layer', 'ml_heat_cap', time_avg=True, files=['atmos_4xdaily']) ## sensible heat flux up at surface, 'W/m^2'

baseexp.use_diag_table(diag_spinup)
baseexp.disable_rrtm()
##baseexp.set_resolution('T85')
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

baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.38 #0.3
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2
baseexp.namelist['spectral_dynamics_nml']['ocean_topog_smoothing'] = 0.8
baseexp.namelist['spectral_dynamics_nml']['raw_filter_coeff'] = 1.
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150.
baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False
baseexp.namelist['spectral_dynamics_nml']['lon_max'] = 256
baseexp.namelist['spectral_dynamics_nml']['lat_max'] = 128
baseexp.namelist['spectral_dynamics_nml']['num_fourier'] = 85
baseexp.namelist['spectral_dynamics_nml']['num_spherical'] = 86


# Land/topography parameters
#heat_cap_ratio = 2./20.
LENKA=[10.0,40.0, True] ## Novak et al 2018 [d,qflux_amp]

baseexp.namelist['mixed_layer_nml']['depth'] = float(depth)
baseexp.namelist['mixed_layer_nml']['load_qflux'] = False
if (qf>0):
    baseexp.namelist['mixed_layer_nml']['do_qflux'] = True
    baseexp.namelist['mixed_layer_nml']['qflux_amp'] = float(qf)
    baseexp.namelist['mixed_layer_nml']['qflux_width'] = 11.3
else :
    baseexp.namelist['mixed_layer_nml']['do_qflux'] = False
baseexp.namelist['mixed_layer_nml']['do_read_sst']=False

## To match with Ruth's parameters

baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray']   = True
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False
baseexp.namelist['idealized_moist_phys_nml']['turb'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_damping'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_virtual'] = False
baseexp.namelist['idealized_moist_phys_nml']['roughness_mom'] = 2e-04
baseexp.namelist['idealized_moist_phys_nml']['roughness_heat'] = 2e-04
baseexp.namelist['idealized_moist_phys_nml']['roughness_moist'] = 2e-04

baseexp.namelist['two_stream_gray_rad_nml']['linear_tau']   = 0.2
baseexp.namelist['two_stream_gray_rad_nml']['wv_exponent']   = 4.0
baseexp.namelist['two_stream_gray_rad_nml']['solar_exponent']   = 2.0
baseexp.namelist['two_stream_gray_rad_nml']['ir_tau_pole']      = 1.8
baseexp.namelist['two_stream_gray_rad_nml']['ir_tau_eq']   = 7.2
baseexp.namelist['two_stream_gray_rad_nml']['atm_abs']     = 0.22
baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal']     = True
baseexp.namelist['two_stream_gray_rad_nml']['dt_rad_avg']      = 86400
baseexp.namelist['two_stream_gray_rad_nml']['use_time_average_coszen'] = True
baseexp.namelist['two_stream_gray_rad_nml']['equinox_day'] = 0.75

## Radiative gases mixixng ratio
baseexp.namelist['constants_nml']['omega'] = 7.2921150e-5

baseexp.runmonth(1, use_restart=False, num_cores=16)
for i in range(2, 361):  
    baseexp.runmonth(i, num_cores=16)
