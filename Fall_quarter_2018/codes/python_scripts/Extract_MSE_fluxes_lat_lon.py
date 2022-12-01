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
log_directory='/project2/tas1/pragallva/Fall_quarter_2017/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

# try:
#    import cPickle as pickle        Unfortunately pickle doesn't work
# except:
#    import pickle

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

#dirc=sys.argv
source='/project2/tas1/pragallva/Fall_quarter_2017/exp_data/'+dirc[1]+'/'
one_year=source+dirc[1]+'_'+dirc[2]+'.nc'

destination='/project2/tas1/pragallva/Fall_quarter_2017/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'
make_sure_path_exists(destination);
v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('...........Extract_ncfile_save_fluxes_radiation.py.............')

logging.debug('imported nc file')

a=6371.0e3  ## m
start=0
end=-1

sigma= 5.670367e-8 ## Stefan's constant
Rd=286.9 # J/Kg
Rv=461.5 # J/kg
Cp= 1004.64 # J/kg/deg
g= 9.8
L=2.500e6 # J/kg
#mon=12

times=v_var['time'][:]
mon=len(times)
sigma_full=v_var['pfull'][::-1]
sigma_half=v_var['phalf'][::-1]
sigma_half[-1]=0.0001
p_sfc=v_var['ps'][:]  ## Pa
lat=v_var['lat'][:]
lon=v_var['lon'][:]

time= nc.num2date(times, units=v_var['time'].units, calendar= v_var['time'].calendar )

logging.debug('imported coordinates')

# In[2]:

pres_full=v_var['pres_full'][:,::-1,:,:]
pres_half=v_var['pres_half'][:,::-1,:,:]

logging.debug('imported pressure')

# In[4]:

time      =nc.num2date(times,units=v_var['time'].units, calendar=v_var['time'].calendar)

# In[10]:

temp      =v_var['temp'][:,::-1,:,:] ## (time, lev, lat, lon)
logging.debug('saved temp')
v_comp    =v_var['vcomp'][:,::-1,:,:]
logging.debug('saved v_comp')
u_comp    =v_var['ucomp'][:,::-1,:,:]
logging.debug('saved u_comp')
Z         =v_var['height'][:,::-1,:,:]
logging.debug('saved Z')
q         =v_var['sphum'][:,::-1,:,:]  ### Specific humidity
logging.debug('saved q')
EKE       =v_var['EKE'][...]  ### Specific humidity
logging.debug('saved EKE')
  
#virtual_T = (1.0+q*(Rv/Rd-1))*temp
#logging.debug('calculated virtual T')
# In[13]:


CpT       =  Cp*temp
Lq        =  L*q
gZ        =  g*Z
# gZ_calc =  g*Z_calc
MSE       = (CpT+Lq+gZ)
MSE_flux  = v_comp*(MSE)

logging.debug('calculated MSE fluxes')


#------------------------------------------------------------

import scipy.integrate as integrate

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(lat))*a
    return y/10**15

def reshape(y):
    y1=y.reshape((12,30,4,len(sigma_full),len(lat),len(lon)))
    return y1

def reshape_pres(y):
    y1=y.reshape((12,30,4,len(sigma_half),len(lat),len(lon)))
    return y1

ps=reshape_pres(pres_half)
weights = (ps[:,:,:,:-1,:,:]-ps[:,:,:,1:,:,:])/g

def weighted(arg):
    w = arg*weights
    return w

def MSE_total(m):
    v=v_comp
    M         = reshape(m)  # (month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    flux      = (V)*weighted(M)            # (month,days, hour, plev, lat, lon)
    monthly_f = flux.mean(axis=1).mean(axis=1)  # (month, plev, lat, lon)
    vert_flux = monthly_f.sum(axis=1)           # (month, lat, lon)
    return vert_flux

def mean_meridional(m):
    v=v_comp
    nlon=len(lon)
    M         = reshape(m)  # (month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = (M).mean(axis=1).mean(axis=1)            # (month,days, hour, plev, lat, lon)
    monthly_v = weighted(V).mean(axis=1).mean(axis=1)  # (month, plev, lat, lon)
    zonal_m   = monthly_m.mean(axis=-1)        # (month, plev, lat)
    zonal_v   = monthly_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=1)           # (month, lat)
    vert_flux = np.dstack([vert_flux]*nlon)
    return vert_flux
    
def SE_vstar_mstar(m):
    v=v_comp
    M         = reshape(m)               # (month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = (M).mean(axis=1).mean(axis=1)           # (month,days, hour, plev, lat, lon)
    monthly_v = (V).mean(axis=1).mean(axis=1)                             # (month, plev, lat, lon)
    m_star    = monthly_m-monthly_m.mean(axis=-1)[...,None]               # (month, plev, lat,lon)
    v_star    = monthly_v-monthly_v.mean(axis=-1)[...,None]
    flux_weighted  = weighted((m_star*v_star)[:,None,None,:,:,:])          # (month, plev, lat, lon)
    vert_flux = flux_weighted.mean(axis=1).mean(axis=1).sum(axis=1)         # (month, lat, lon)
    return vert_flux

def SE_vbar_mstar(m):
    v=v_comp
    M         = reshape(m)               # (month,days, hour, plev, lat, lon)
    V         = reshape(v)               # (month,days, hour, plev, lat, lon)
    monthly_Zv =(V).mean(axis=1).mean(axis=1).mean(axis=-1)                             # (month, plev, lat)
    monthly_m = (M).mean(axis=1).mean(axis=1)
    m_star    = monthly_m-monthly_m.mean(axis=-1)[...,None]               # (month, plev, lat,lon)
    flux_weighted  = weighted((monthly_Zv[...,None]*m_star)[:,None,None,:,:,:])          # (month, plev, lat, lon)
    vert_flux = flux_weighted.mean(axis=1).mean(axis=1).sum(axis=1)         # (month, lat, lon)
    return vert_flux

def SE_vstar_mbar(m):
    v=v_comp
    M         = reshape(m)               # (month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_Zm = (M).mean(axis=1).mean(axis=1).mean(axis=-1)               # (month, plev, lat)
    monthly_v = (V).mean(axis=1).mean(axis=1)              # (month, plev, lat)
    m_star    = monthly_m-monthly_m.mean(axis=-1)[...,None]               # (month, plev, lat,lon)
    v_star    = monthly_v-monthly_v.mean(axis=-1)[...,None]
    flux_weighted  = weighted((monthly_Zm[...,None]*v_star)[:,None,None,:,:,:])          # (month, plev, lat, lon)
    vert_flux = flux_weighted.mean(axis=1).mean(axis=1).sum(axis=1)         # (month, lat, lon)
    return vert_flux
    
def transient_eddies(m):
    v=v_comp
    M         = reshape(m)                       # (month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = M.mean(axis=1).mean(axis=1)                     # average over days 
    monthly_v = V.mean(axis=1).mean(axis=1)           # (month, plev, lat, lon)
    m_prime   = M-monthly_m[:,None,None,:,:,:]                  # (month, days, plev, lat, lon)
    v_prime   = V-monthly_v[:,None,None,:,:,:]
    flux      = weighted(m_prime*v_prime).mean(axis=1).mean(axis=1)  # (month,  plev, lat, lon) 
    vert_flux = (flux).sum(axis=1)                        # (month,  lat, lon) 
    return vert_flux

logging.debug('Defined functions')

MSE_flux=MSE_total(CpT+gZ+Lq)
MM_flux=mean_meridional(CpT+gZ+Lq)
SE1_flux=SE_vstar_mstar(CpT+gZ+Lq)
SE2_flux=SE_vbar_mstar(CpT+gZ+Lq)
SE3_flux=SE_vstar_mbar(CpT+gZ+Lq)
TE_flux=transient_eddies(CpT+gZ+Lq)


fluxes_dic      ={"MSE_flux":MSE_flux,"MM_flux":MM_flux,'SE1_flux':SE1_flux,'SE2_flux':SE2_flux,'SE3_flux':SE3_flux}


logging.debug("loaded lat lon MSE_flux dictionary")

# SAVING AS A DICTIONARY

plev = np.array([0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0])

coord_dic       ={"lat":lat,"lon":lon,"time":times,"plev": plev} 

logging.debug("loaded coordinates dictionary")


# In[141]:
save(destination+"coord_dic.hkl"       ,coord_dic)
save(destination+"fluxes_latlon_dic.hkl",fluxes_dic)
#save(destination+"rad_dic.hkl"         ,rad_dic)

logging.debug("successfully saved lat lon dictionaries in files")


