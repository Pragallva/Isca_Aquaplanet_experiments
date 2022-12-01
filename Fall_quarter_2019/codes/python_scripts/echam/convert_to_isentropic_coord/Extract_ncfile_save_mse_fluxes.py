# coding: utf-8

import netCDF4 as nc
import pylab as py
import numpy as np
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.cm as cm
import sys
import os
import errno
dirc=sys.argv

import logging
import time as ti
start = ti.time()

CO2         = str(dirc[1])
num         = str(dirc[2])
no_of_days  = 30

current_dir = os.getcwd()
parts       = current_dir.split(os.sep)
quarter     = "/"+os.path.join(*parts[:5])

import hickle as hkl

## A function to save a dictionary ##
def save(filename,dictionary):
    hkl.dump(dictionary, filename, mode='w')
    
## A function to load a dictionary ## 
def load(filename):
    dictionary = hkl.load(filename)
    return dictionary

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
        logging.debug('destination folder created !')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logging.debug('destination folder exists already!')
            raise

log_directory = quarter+'/codes/shell_scripts/convert_isentrope/log/'
make_sure_path_exists( path   = log_directory )


for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

if CO2=='None':
     CO2 = 'echr0001'
exp_dir = 'ATM_dm_pl_'+CO2+'_10'

source = quarter+'/exp_data/echam/data_in_isentropic_coord/'+exp_dir+'/'
make_sure_path_exists(source);
one_year = source+exp_dir+num+'.nc'

logging.basicConfig( filename = log_directory+'momentum'+exp_dir+num+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

destination = quarter+'/post_process_data/echam/data_in_isentropic_coord/'+exp_dir+num+'/'
make_sure_path_exists(destination);

v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('....Extract ncfile save MSE and momentum fluxes and radiation....'+str(num)+' averaged over '+str(no_of_days))

logging.debug('imported nc file')

a=6371.0e3  ## m

sigma= 5.670367e-8 ## Stefan's constant
Rd=286.9 # J/Kg
Rv=461.5 # J/kg
Cp= 1004.64 # J/kg/deg
G= 9.8
L=2.500e6 # J/kg
#mon=12

last = 360

lat=v_var['lat'][:]
lon=v_var['lon'][:]
theta=v_var['level'][:]


logging.debug('imported coordinates')


if no_of_days == None:
   no_of_days = 30
else:
   days = int(no_of_days)

no_of_years = 360/(12*30)
no_of_months =12
hours = 1


pres      = v_var['pres'][:last,...]*100
logging.debug('imported pressure')

temp      = v_var['t'][:last,...] ## (time, lev, lat, lon)
logging.debug('saved temp')

v_comp    = v_var['v'][:last,...]
logging.debug('saved v_comp')

u_comp    = v_var['u'][:last,...]
logging.debug('saved u_comp')

omega     = v_var['omega'][:last,...]
logging.debug('saved omega')

Z         = v_var['geopoth'][:last,...]
logging.debug('saved Z')

q         = v_var['q'][:last,...]  ### Specific humidity
logging.debug('saved q')

theta_dot = v_var['theta_dot'][:last,...]
logging.debug('saved theta_dot')

dtheta_dp = v_var['dtheta_dp'][:last,...]
logging.debug('saved dtheta_dp')

deltheta_delt = v_var['deltheta_delt'][:last,...]
logging.debug('saved deltheta_delt')

rho_s = -(1/G)*( np.gradient(pres,theta,axis=1) )

dgZ_dtheta = np.gradient(G*Z,theta,axis=1)

u_sq= u_comp**2
v_sq= v_comp**2
KE  = (u_sq+v_sq)/2.0

logging.debug('calculated KE')

CpT       =  Cp*temp
Lq        =  L*q
gZ        =  G*Z
MSE       = (CpT+Lq+gZ+KE)

logging.debug('calculated MSE fluxes')


def reshape(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(theta),len(lat),len(lon)))
    return y1

def weighted(arg):
    w = arg*reshape(rho_s)
    return w

def R(y):
    y1 = reshape(y)
    return np.nanmean( np.nanmean( np.nanmean( np.nanmean(reshape(y1), axis=2), axis=2), axis=-1), axis=0) 

logging.debug("loaded coordinates dictionary")

import scipy.integrate as integrate

def MSE_total(m,vv=v_comp):
    v=vv
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)
    flux      = weighted(V)*(M)                # (year, month,days, hour, plev, lat, lon)
    monthly_f = np.nanmean(np.nanmean(flux,axis=2),axis=2) # (year, month, plev, lat, lon)
    zonal_f   = np.nanmean(monthly_f,axis=-1)        # (year, month, plev, lat)
    return np.nanmean(zonal_f,axis=0)


def mean_meridional(m,vv=v_comp):
    v         = vv
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)
    monthly_m = np.nanmean(np.nanmean(M,axis=2),axis=2)  # (year, month, days, hour, plev, lat, lon)
    monthly_v = np.nanmean(np.nanmean(weighted(V),axis=2),axis=2)  # (year, month, plev, lat, lon)
    zonal_m   = np.nanmean(monthly_m,axis=-1)        # (year, month, plev, lat)
    zonal_v   = np.nanmean(monthly_v,axis=-1)
    return np.nanmean(zonal_m*zonal_v, axis=0)


def stationary_eddies(m,vv=v_comp):
    v         = vv
    M         = reshape(m)                              # (year, month, days, hour, plev, lat, lon)
    V         = reshape(v)
    monthly_m = np.nanmean(np.nanmean(M,axis=2),axis=2)           # (year, month,days, hour, plev, lat, lon)
    monthly_v = np.nanmean(np.nanmean(weighted(V),axis=2),axis=2)                        # (year, month, plev, lat, lon)
    m_star    = monthly_m-np.nanmean(monthly_m,axis=-1)[...,None]          # (year, month, plev, lat, lon)
    v_star    = monthly_v-np.nanmean(monthly_v,axis=-1)[...,None]
    flux_weighted  = ((m_star*v_star)[:,:,None,None,:,:,:])    # (year, month, plev, lat)
    return np.nanmean(np.nanmean(np.nanmean(np.nanmean(flux_weighted,axis=2),axis=2),axis=-1),axis=0)


def transient_eddies(m,vv=v_comp):
    v         = vv
    M         = reshape(m)                            # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)
    monthly_m = np.nanmean(np.nanmean(M,axis=2),axis=2)           # average over days 
    monthly_v = np.nanmean(np.nanmean(weighted(V),axis=2),axis=2)           # (month, plev, lat, lon)
    m_prime   = M-monthly_m[:,:,None,None,:,:,:]      # (month, days, plev, lat, lon)
    v_prime   = V-monthly_v[:,:,None,None,:,:,:]
    flux      = np.nanmean(np.nanmean(np.nanmean(m_prime*v_prime,axis=2),axis=2),axis=-1) # (year,month,plev,lat) 
    return np.nanmean(flux,axis=0)


logging.debug("defined functions")

#############################
########   dm_by_dt    ######
#############################

def tendency_24_hours(h):
    dt=24*60*60
    d_by_dt=np.gradient(h,dt,axis=0)
    return d_by_dt

moist_enthalpy  = CpT+Lq+KE   ## ((1440, 29, 64, 128))
dh_by_dt    = R(tendency_24_hours(rho_s*moist_enthalpy))
dm_by_dt    = R(tendency_24_hours(rho_s*MSE))

dgZ_by_dt_term_m1  = R(rho_s*tendency_24_hours(gZ))
dgZ_by_dt_term2    = R(rho_s*deltheta_delt*dgZ_dtheta) 

dgZ_by_dt_term_h1  = -R(gZ*tendency_24_hours(rho_s))


logging.debug("calculated tendency dictionary")

mse_flux  = MSE_total(CpT+gZ+Lq+KE)

######## divergence function ######

def div_y(X):
   cos_phi = np.cos(np.deg2rad(lat))[None,None,:]
   div     = np.gradient(X*cos_phi,np.deg2rad(lat),axis=-1)/(a*cos_phi)
   return div

def div_w(X):
   div    = np.gradient(X, theta, axis=1)
   return div 

logging.debug("defined divergence functions")
######## v direction #######

mse_flux_v  = MSE_total(MSE) 

MM_flux_v     =mean_meridional(MSE)
SE_flux_v     =stationary_eddies(MSE)
TE_flux_v     =transient_eddies(MSE)

logging.debug("Total fluxes v")

TE_sensible_v=transient_eddies(CpT)
TE_moist_v   =transient_eddies(Lq)
TE_pot_v     =transient_eddies(gZ)
TE_KE_v      =transient_eddies(KE)

logging.debug("TE fluxes v")

SE_sensible_v=stationary_eddies(CpT)
SE_moist_v   =stationary_eddies(Lq)
SE_pot_v     =stationary_eddies(gZ)
SE_KE_v      =stationary_eddies(KE)

logging.debug("SE fluxes v")

MM_sensible_v=mean_meridional(CpT)
MM_moist_v   =mean_meridional(Lq)
MM_pot_v     =mean_meridional(gZ)
MM_KE_v      =mean_meridional(KE)

logging.debug("MM fluxes v")

v_fluxes= {'MM_flux_v':MM_flux_v,         'SE_flux_v':SE_flux_v,   'TE_flux_v':TE_flux_v, 'mse_flux_v': mse_flux_v,\
           'TE_sensible_v':TE_sensible_v, 'TE_moist_v':TE_moist_v, 'TE_pot_v':TE_pot_v, 'TE_KE_v':TE_KE_v,\
           'SE_sensible_v':SE_sensible_v, 'SE_moist_v':SE_moist_v, 'SE_pot_v':SE_pot_v, 'SE_KE_v':SE_KE_v,\
           'MM_sensible_v':MM_sensible_v, 'MM_moist_v':MM_moist_v, 'MM_pot_v':MM_pot_v, 'MM_KE_v':MM_KE_v}

logging.debug("v flux dictionary")

######## w direction #######
mse_flux_w  = MSE_total(MSE, theta_dot)

MM_flux_w     =mean_meridional(MSE,theta_dot)
SE_flux_w     =stationary_eddies(MSE,theta_dot)
TE_flux_w     =transient_eddies(MSE,theta_dot)

logging.debug("Total fluxes w")

TE_sensible_w=transient_eddies(CpT,theta_dot)
TE_moist_w   =transient_eddies(Lq,theta_dot)
TE_pot_w     =transient_eddies(gZ,theta_dot)
TE_KE_w      =transient_eddies(KE,theta_dot)

logging.debug("TE fluxes w")

SE_sensible_w=stationary_eddies(CpT,theta_dot)
SE_moist_w   =stationary_eddies(Lq,theta_dot)
SE_pot_w     =stationary_eddies(gZ,theta_dot)
SE_KE_w      =stationary_eddies(KE,theta_dot)


logging.debug("SE fluxes w")

MM_sensible_w=mean_meridional(CpT,theta_dot)
MM_moist_w   =mean_meridional(Lq,theta_dot)
MM_pot_w     =mean_meridional(gZ,theta_dot)
MM_KE_w      =mean_meridional(KE,theta_dot)

logging.debug("MM fluxes w")

w_fluxes= {'MM_flux_w':MM_flux_w,         'SE_flux_w':SE_flux_w,   'TE_flux_w':TE_flux_v, 'mse_flux_w':mse_flux_w,\
           'TE_sensible_w':TE_sensible_w, 'TE_moist_w':TE_moist_w, 'TE_pot_w':TE_pot_v, 'TE_KE_w':TE_KE_v,\
           'SE_sensible_w':SE_sensible_w, 'SE_moist_w':SE_moist_w, 'SE_pot_w':SE_pot_v, 'SE_KE_w':SE_KE_v,\
           'MM_sensible_w':MM_sensible_w, 'MM_moist_w':MM_moist_w, 'MM_pot_w':MM_pot_v, 'MM_KE_w':MM_KE_v}

logging.debug("w flux dictionary")
##############################


######## div_y direction #######

div_mse_flux_v    = div_y(mse_flux_v)

div_MM_flux_v     = div_y(MM_flux_v)
div_SE_flux_v     = div_y(SE_flux_v)
div_TE_flux_v     = div_y(TE_flux_v)

div_TE_sensible_v= div_y(TE_sensible_v)
div_TE_moist_v   = div_y(TE_moist_v)
div_TE_pot_v     = div_y(TE_pot_v)
div_TE_KE_v      = div_y(TE_KE_v)

div_SE_sensible_v= div_y(SE_sensible_v)
div_SE_moist_v   = div_y(SE_moist_v)
div_SE_pot_v     = div_y(SE_pot_v)
div_SE_KE_v      = div_y(SE_KE_v)

div_MM_sensible_v= div_y(MM_sensible_v)
div_MM_moist_v   = div_y(MM_moist_v)
div_MM_pot_v     = div_y(MM_pot_v)
div_MM_KE_v      = div_y(MM_KE_v)

logging.debug("Div v")

div_v_fluxes= {'div_MM_flux_v':div_MM_flux_v, 'div_SE_flux_v':div_SE_flux_v,   'div_TE_flux_v':div_TE_flux_v, 'div_mse_flux_v':div_mse_flux_v,\
         'div_TE_sensible_v':div_TE_sensible_v, 'div_TE_moist_v':div_TE_moist_v, 'div_TE_pot_v':div_TE_pot_v, 'div_TE_KE_v':div_TE_KE_v,\
         'div_SE_sensible_v':div_SE_sensible_v, 'div_SE_moist_v':div_SE_moist_v, 'div_SE_pot_v':div_SE_pot_v, 'div_SE_KE_v':div_SE_KE_v,\
         'div_MM_sensible_v':div_MM_sensible_v, 'div_MM_moist_v':div_MM_moist_v, 'div_MM_pot_v':div_MM_pot_v, 'div_MM_KE_v':div_MM_KE_v}

logging.debug("Div v dictionary")

######## div_w direction #######

div_mse_flux_w    = div_y(mse_flux_w)

div_MM_flux_w     = div_w(MM_flux_w)
div_SE_flux_w     = div_w(SE_flux_w)
div_TE_flux_w     = div_w(TE_flux_w)

div_TE_sensible_w= div_w(TE_sensible_w)
div_TE_moist_w   = div_w(TE_moist_w)
div_TE_pot_w     = div_w(TE_pot_w)
div_TE_KE_w      = div_w(TE_KE_w)

div_SE_sensible_w= div_w(SE_sensible_w)
div_SE_moist_w   = div_w(SE_moist_w)
div_SE_pot_w     = div_w(SE_pot_w)
div_SE_KE_w      = div_w(SE_KE_w)

div_MM_sensible_w= div_w(MM_sensible_w)
div_MM_moist_w   = div_w(MM_moist_w)
div_MM_pot_w     = div_w(MM_pot_w)
div_MM_KE_w      = div_w(MM_KE_w)

logging.debug("Div w")

div_w_fluxes= {'div_MM_flux_w':div_MM_flux_w, 'div_SE_flux_w':div_SE_flux_w, 'div_TE_flux_w':div_TE_flux_w,  'div_mse_flux_w':div_mse_flux_w,\
         'div_TE_sensible_w':div_TE_sensible_w, 'div_TE_moist_w':div_TE_moist_w, 'div_TE_pot_w':div_TE_pot_w, 'div_TE_KE_w':div_TE_KE_w,\
         'div_SE_sensible_w':div_SE_sensible_w, 'div_SE_moist_w':div_SE_moist_w, 'div_SE_pot_w':div_SE_pot_w, 'div_SE_KE_w':div_SE_KE_w,\
         'div_MM_sensible_w':div_MM_sensible_w, 'div_MM_moist_w':div_MM_moist_w, 'div_MM_pot_w':div_MM_pot_w, 'div_MM_KE_w':div_MM_KE_w}

logging.debug("Div w dictionary")
#################################

logging.debug("calculated tendencies")

fluxes_PW  = {'v_fluxes':v_fluxes, 'w_fluxes':w_fluxes}
fluxes_Wm2 = {'div_v_fluxes':div_v_fluxes, 'div_w_fluxes':div_w_fluxes}

tendency_terms = {'dh_by_dt':dh_by_dt, 'dm_by_dt':dm_by_dt, 'dgZ_by_dt_term_m1':dgZ_by_dt_term_m1,\
                  'dgZ_by_dt_term2':dgZ_by_dt_term2,        'dgZ_by_dt_term_h1':dgZ_by_dt_term_h1}

raw_data      = {'MSE':R(MSE), 'CpT':R(CpT), 'gZ':R(gZ), 'Lq':R(Lq), 'rho_s':R(rho_s),\
                 'u':R(u_comp), 'v':R(v_comp), 'omega':R(omega), 'pres':R(pres)}

coord        = {'lat':lat, 'lon': lon, 'theta': theta}


save(destination+"fluxes_PW.hkl",fluxes_PW)
logging.debug("loaded fluxes_PW dictionary")

save(destination+"fluxes_Wm2.hkl",fluxes_Wm2)
logging.debug("loaded fluxes_Wm2 dictionary")

save(destination+"tendency_terms.hkl",tendency_terms)
logging.debug("loaded tendency terms dictionary")

save(destination+"raw_data.hkl",raw_data)
logging.debug("loaded raw data dictionary")

save(destination+"coord.hkl",coord)
logging.debug("loaded coord data dictionary")

end = ti.time()
logging.debug("Time taken --> "+ str(end-start))
logging.debug("Awesome! complete!")
logging.debug("-------------------------------")




