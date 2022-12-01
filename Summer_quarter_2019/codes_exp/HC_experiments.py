import numpy as np
import os
import math as mth
import sys
dirc=sys.argv

from gfdl import Experiment, DiagTable
from gfdl import land_generator_fn as lgf

edge        = int(dirc[1]) #tropical edge
land_depth  = int(dirc[2])
ocean_depth = int(dirc[3])

experiment_name = 'HC'+str(edge)+'_la'+str(land_depth)+'m'+'_oc'+str(ocean_depth)+'m' #landice55_z0
land_name       = 'HC'+str(edge)

lgf.write_land('HC'+str(edge),land_mode='square', boundaries=[-1*float(edge),float(edge),0.,360.], continents=None, topo_mode= None, mountains=None, topo_gauss=None, waterworld=True)

baseexp         = Experiment(experiment_name, overwrite_data=False)

baseexp.inputfiles = ['/project/tas1/pragallva/Fall_quarter_2017/Isca/input/ozone_1990.nc','/project2/tas1/pragallva/Summer_quarter_2018/land_files/HC'+str(edge)+'.nc']
   
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

diag_spinup.add_field('atmosphere', 'precipitation',  time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('atmosphere', 'tau_u',          time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('atmosphere', 'tau_v',          time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('atmosphere', 'drag_m',         time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('atmosphere', 'w_atm',          time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('atmosphere', 'ustar',          time_avg=True, files=['atmos_4xdaily'])

diag_spinup.add_field('damping',    'udt_rdamp',      time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('damping',    'vdt_rdamp',      time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('damping',    'tdt_diss_rdamp', time_avg=True, files=['atmos_4xdaily'])
diag_spinup.add_field('damping',    'diss_heat_rdamp',time_avg=True, files=['atmos_4xdaily'])

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


baseexp.use_diag_table(diag_spinup)

baseexp.compile()

baseexp.clear_rundir()


baseexp.namelist['main_nml'] = {
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
}

baseexp.namelist['mixed_layer_nml']['depth'] = ocean_depth # float(dirc[1])
baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.25 #0.3
baseexp.namelist['mixed_layer_nml']['tconst'] = 285.
baseexp.namelist['spectral_dynamics_nml']['surf_res'] = 0.2
baseexp.namelist['spectral_dynamics_nml']['ocean_topog_smoothing'] = 0.8
baseexp.namelist['spectral_dynamics_nml']['raw_filter_coeff'] = 1.
baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 150.
baseexp.namelist['rrtm_radiation_nml']['dt_rad'] = 3600
baseexp.namelist['rrtm_radiation_nml']['do_read_ozone']=True
baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False

# Land/topography parameters

heat_cap_ratio=float(land_depth)/float(ocean_depth)
baseexp.namelist['mixed_layer_nml']['land_h_capacity_prefactor'] = heat_cap_ratio
baseexp.namelist['mixed_layer_nml']['land_albedo_prefactor']     = 1 #land_albedo_prefactor # 1.3 #1.2
baseexp.namelist['surface_flux_nml']['land_humidity_prefactor']  = 1 #land_humidity_prefactor # 1.3 #1.2
baseexp.namelist['mixed_layer_nml']['load_qflux'] = False
baseexp.namelist['mixed_layer_nml']['do_qflux'] = False
baseexp.namelist['mixed_layer_nml']['update_albedo_from_ice'] = False

## All these land options are not required 
baseexp.namelist['idealized_moist_phys_nml']['land_option'] = 'input'
baseexp.namelist['mixed_layer_nml']['land_option_zonal']    = False
baseexp.namelist['idealized_moist_phys_nml']['land_file_name'] = 'INPUT/'+land_name+'.nc'
baseexp.namelist['mixed_layer_nml']['land_option'] = 'input'
baseexp.namelist['spectral_init_cond_nml']['topography_option'] = 'flat'#'input'
#baseexp.namelist['spectral_init_cond_nml']['topog_file_name'] = land_name+'.nc'

## To match with Ruth's parameters

baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = False
baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = True
baseexp.namelist['idealized_moist_phys_nml']['turb'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_damping'] = True
baseexp.namelist['idealized_moist_phys_nml']['mixed_layer_bc'] = True
baseexp.namelist['idealized_moist_phys_nml']['do_virtual'] = False
baseexp.namelist['idealized_moist_phys_nml']['roughness_mom'] = 2e-04
baseexp.namelist['idealized_moist_phys_nml']['roughness_heat'] = 2e-04
baseexp.namelist['idealized_moist_phys_nml']['roughness_moist'] = 2e-04

## Radiative gases mixixng ratio
baseexp.namelist['rrtm_radiation_nml']['include_secondary_gases'] = True
baseexp.namelist['rrtm_radiation_nml']['co2ppmv'] = 348.0
baseexp.namelist['rrtm_radiation_nml']['ch4_val'] = 1650.e-9
baseexp.namelist['rrtm_radiation_nml']['n2o_val'] = 306e-9
baseexp.namelist['rrtm_radiation_nml']['cfc11_val'] = 124.7e-12
baseexp.namelist['rrtm_radiation_nml']['cfc12_val'] = 232.9e-12

baseexp.runmonth(1, use_restart=False, num_cores=16)
for i in range(2, 481):  
    baseexp.runmonth(i, num_cores=16)
