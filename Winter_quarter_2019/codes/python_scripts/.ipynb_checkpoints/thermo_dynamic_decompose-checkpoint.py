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
# aqua isca5m

edge =str(dirc[1])
land =str(dirc[2])
ocean=str(dirc[3]) 
num  =str(dirc[4]) ## represents an annual year of data

log_directory='/project2/tas1/pragallva/Winter_quarter_2019/codes/shell_script/log/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m_MMC_decomp'
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

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


source='/project2/tas1/pragallva/Winter_quarter_2019/exp_data/precip/HC'+edge+'_la'+land+'m_oc'+ocean+'m/'
one_year=source+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+num+'.nc'

destination='/project2/tas1/pragallva/Winter_quarter_2019/post_process_data/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+num+'/'
make_sure_path_exists(destination);
v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('...........Extract_ncfile_save_fluxes_radiation.py.............'+str(num))

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


no_of_years = len(time)/(12*30*4)
no_of_months =12
days = 30
hours = 4

# In[2]:

pres_full=v_var['pres_full'][:,::-1,:,:]
pres_half=v_var['pres_half'][:,::-1,:,:]

logging.debug('imported pressure')

# In[4]:

#time      =nc.num2date(times,units=v_var['time'].units, calendar=v_var['time'].calendar)

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



#logging.debug('calculated KE')

CpT       =  Cp*temp
Lq        =  L*q
gZ        =  g*Z
# gZ_calc =  g*Z_calc

logging.debug('calculated MSE fluxes')

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)


logging.debug('saved raw data')



# In[138]:

################################################################ 
logging.debug("Begins calculating TE, SE and MM fluxes ....")

import scipy.integrate as integrate
def integrated(x):
    l=np.deg2rad(lat)
    x=x*np.cos(l)
    int_x  =integrate.cumtrapz(x[::-1],l[::-1],axis=0,initial=None) #  (This is basically integration from - 90 deg)
    int_x_r=integrate.cumtrapz(x      ,l      ,axis=0,initial=None) #  (This is basically integration from + 90 deg) 
    avg_int_r=2*np.pi*a**2*(int_x[::-1][1:]+int_x_r[:-1])/2.0
    return avg_int_r/10**15

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(lat))*a
    return y/10**15

def reshape(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1

def reshape_pres(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_half),len(lat),len(lon)))
    return y1

ps=reshape_pres(pres_half)
weights = (ps[:,:,:,:,:-1,:,:]-ps[:,:,:,:,1:,:,:])/g

def weighted(arg):
    w = arg*weights
    return w
    
def mean_meridional(m,vert=1):
    v=v_comp
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = (M).mean(axis=2).mean(axis=2)  # (year, month, days, hour, plev, lat, lon)
    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    zonal_m   = monthly_m.mean(axis=-1)        # (year, month, plev, lat)
    zonal_v   = monthly_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=2)  # (year, month, lat)
    if vert==1:
        return vert_flux.mean(axis=0)              # (month, lat)--> yearly average
    else :
        return (zonal_m*zonal_v).mean(axis=0)
    
def delta_v_mmc(m,control,control_vert,vert=1):
    v         = v_comp
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)

    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    monthly_m = (M).mean(axis=2).mean(axis=2)

    delta_v   = monthly_v-monthly_v.mean(axis=1)[:,None,...]   # (year, month, plev, lat, lon)
    average_m = monthly_m.mean(axis=1)[:,None,...] # (year, month, plev, lat, lon)

    zonal_m   = average_m.mean(axis=-1)                        # (year, month, plev, lat)
    zonal_v   = delta_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=2)  # (year, month, lat)
    if vert==1:
        return (vert_flux).mean(axis=0)# (month, lat)--> yearly average
    else :
        return (zonal_m*zonal_v).mean(axis=0)
    

def delta_m_mmc(m,control,control_vert,vert=1):
    v         = v_comp
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)

    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    monthly_m = (M).mean(axis=2).mean(axis=2)

    delta_m   = monthly_m-monthly_m.mean(axis=1)[:,None,...]   # (year, month, plev, lat, lon)
    average_v = monthly_v.mean(axis=1)[:,None,...] # (year, month, plev, lat,lon)

    zonal_m   = delta_m.mean(axis=-1)                        # (year, month, plev, lat)
    zonal_v   = average_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=2)  # (year, month, lat)
    if vert==1:
        return (vert_flux).mean(axis=0)# (month, lat)--> yearly average
    else :
        return (zonal_m*zonal_v).mean(axis=0)
    

def delta_mv_mmc(m,control,control_vert,vert=1):
    v         = v_comp
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)

    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    monthly_m = (M).mean(axis=2).mean(axis=2)

    delta_m   = monthly_m-monthly_m.mean(axis=1)[:,None,...]   # (year, month, plev, lat, lon)
    delta_v   = monthly_v-monthly_v.mean(axis=1)[:,None,...] # (year, month, plev, lat,lon)

    zonal_m   = delta_m.mean(axis=-1)                        # (year, month, plev, lat)
    zonal_v   = delta_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=2)  # (year, month, lat)
    if vert==1:
        return (vert_flux).mean(axis=0)# (month, lat)--> yearly average
    else :
        return (zonal_m*zonal_v).mean(axis=0)

#moist_enthalpy = CpT+Lq+KE   ## ((1440, 29, 64, 128))

MM_flux    =mean_meridional(CpT+gZ+Lq)
#MM_sensible=mean_meridional(CpT)
#MM_moist   =mean_meridional(Lq)
#MM_pot     =mean_meridional(gZ)

del_m_mmc           =delta_m_mmc(CpT+gZ+Lq)
del_m_mmc_sensible  =delta_m_mmc(CpT)
del_m_mmc_moist     =delta_m_mmc(Lq)
del_m_mmc_pot       =delta_m_mmc(gZ)

del_v_mmc           =delta_v_mmc(CpT+gZ+Lq )
del_v_mmc_sensible  =delta_v_mmc(CpT )
del_v_mmc_moist     =delta_v_mmc(Lq )
del_v_mmc_pot       =delta_v_mmc(gZ )

del_mv_mmc           =delta_mv_mmc(CpT+gZ+Lq )
del_mv_mmc_sensible  =delta_mv_mmc(CpT )
del_mv_mmc_moist     =delta_mv_mmc(Lq )
del_mv_mmc_pot       =delta_mv_mmc(gZ )

##### Few more linesof code

del_m_mmc_vert           =delta_m_mmc(CpT+gZ+Lq ,vert=0)
del_m_mmc_sensible_vert  =delta_m_mmc(CpT,vert=0)
del_m_mmc_moist_vert     =delta_m_mmc(Lq, vert=0)
del_m_mmc_pot_vert       =delta_m_mmc(gZ, vert=0)

del_v_mmc_vert           =delta_v_mmc(CpT+gZ+Lq,vert=0)
del_v_mmc_sensible_vert  =delta_v_mmc(CpT, vert=0)
del_v_mmc_moist_vert     =delta_v_mmc(Lq, vert=0)
del_v_mmc_pot_vert       =delta_v_mmc(gZ,vert=0)

del_mv_mmc_vert           =delta_mv_mmc(CpT+gZ+Lq,vert=0)
del_mv_mmc_sensible_vert  =delta_mv_mmc(CpT, vert=0)
del_mv_mmc_moist_vert     =delta_mv_mmc(Lq, vert=0)
del_mv_mmc_pot_vert       =delta_mv_mmc(gZ,vert=0)


MMC_thermo_dynamic_dic   ={"MM_flux":MM_flux,"del_m_mmc":del_m_mmc, "del_v_mmc":del_v_mmc , "del_m_mmc_sensible":del_m_mmc_sensible, "del_v_mmc_sensible":del_v_mmc_sensible, "del_v_mmc_moist":del_v_mmc_moist,"del_m_mmc_moist":del_m_mmc_moist,"del_v_mmc_pot":del_v_mmc_pot,"del_m_mmc_pot":del_m_mmc_pot,"del_mv_mmc":del_mv_mmc , "del_mv_mmc_sensible":del_mv_mmc_sensible, "del_mv_mmc_moist":del_mv_mmc_moist, "del_mv_mmc_pot":del_mv_mmc_pot}

save(destination+"MMC_thermo_dynamic_dic.hkl" ,MMC_thermo_dynamic_dic)
logging.debug("loaded MMC_thermo_dynamic dictionary")

MMC_thermo_dynamic_dic_vert   ={"MM_flux":MM_flux_vert,"del_m_mmc":del_m_mmc_vert, "del_v_mmc":del_v_mmc_vert , "del_m_mmc_sensible":del_m_mmc_sensible_vert, "del_v_mmc_sensible":del_v_mmc_sensible_vert, "del_v_mmc_moist":del_v_mmc_moist_vert,"del_m_mmc_moist":del_m_mmc_moist_vert,"del_v_mmc_pot":del_v_mmc_pot_vert,"del_m_mmc_pot":del_m_mmc_pot_vert,"del_mv_mmc":del_mv_mmc_vert , "del_mv_mmc_sensible":del_mv_mmc_sensible_vert, "del_mv_mmc_moist":del_mv_mmc_moist_vert, "del_mv_mmc_pot":del_mv_mmc_pot_vert}

save(destination+"MMC_thermo_dynamic_dic_vert.hkl" ,MMC_thermo_dynamic_dic_vert)
logging.debug("loaded MMC_thermo_dynamic_vert dictionary")

