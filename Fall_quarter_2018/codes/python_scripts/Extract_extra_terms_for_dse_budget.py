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

edge=str(dirc[1])
land=str(dirc[2])
ocean=str(dirc[3])

num=str(dirc[4]) ## represents an annual year of data

log_directory='/project2/tas1/pragallva/Fall_quarter_2018/codes/shell_script/log/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'
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


source='/project2/tas1/pragallva/Fall_quarter_2018/exp_data/precip/HC'+edge+'_la'+land+'m_oc'+ocean+'m/'
one_year=source+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+num+'.nc'

destination='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+num+'/'
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

pres_full=v_var['pres_full'][:,::-1,:,:]
pres_half=v_var['pres_half'][:,::-1,:,:]

logging.debug('imported pressure')

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

u_sq= u_comp**2
v_sq= v_comp**2
KE  = (u_sq+v_sq)/2.0

logging.debug('calculated KE')

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)

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

def tendency_six_hours(h):
    dt=6*60*60
    dh_by_dt=np.copy(h)
    for lo in range(len(lon)):
        for la in range(len(lat)):
            for lev in range(len(sigma_full)):
                dh_by_dt[:,lev,la,lo]=np.gradient( h[:,lev,la,lo],dt)
    return dh_by_dt

logging.debug("calculated dh_by_dt dictionary")

CpT = Cp*temp
Lq  = L*q
gZ  = g*Z

moist_enthalpy = CpT+Lq+KE   ## ((1440, 29, 64, 128))
dh_by_dt=tendency_six_hours(moist_enthalpy)
dsens_by_dt=tendency_six_hours(CpT)
dpot_by_dt=tendency_six_hours(gZ)
dmoist_by_dt=tendency_six_hours(Lq)


dhdt=weighted(reshape(dh_by_dt)).mean(axis=2).mean(axis=2).mean(axis=-1) # reshape 
dhdt_vert=dhdt.sum(axis=2).mean(axis=0)  # (month, lat)--> yearly average and vertical average
dhdt_vert_vert=dhdt.mean(axis=0)  ## No vertical integration

dsensdt=weighted(reshape(dsens_by_dt)).mean(axis=2).mean(axis=2).mean(axis=-1) # reshape 
dsensdt_vert=dsensdt.sum(axis=2).mean(axis=0)  # (month, lat)--> yearly average and vertical average
dsensdt_vert_vert=dsensdt.mean(axis=0)  ## No vertical integration

dmoistdt=weighted(reshape(dmoist_by_dt)).mean(axis=2).mean(axis=2).mean(axis=-1) # reshape 
dmoistdt_vert=dmoistdt.sum(axis=2).mean(axis=0)  # (month, lat)--> yearly average and vertical average
dmoistdt_vert_vert=dmoistdt.mean(axis=0)  ## No vertical integration

dpotdt=weighted(reshape(dpot_by_dt)).mean(axis=2).mean(axis=2).mean(axis=-1) # reshape 
dpotdt_vert=dpotdt.sum(axis=2).mean(axis=0)  # (month, lat)--> yearly average and vertical average
dpotdt_vert_vert=dpotdt.mean(axis=0)  ## No vertical integration

logging.debug("calculated dh_by_dt and individual tendency terms for sensible, moist and potential")

# SAVING AS A DICTIONARY
tendency_Wm2_dic = {"dhdt":dhdt_vert,"dhdt_vert":dhdt_vert,'dsensdt':dsensdt_vert,'dsensdt_vert':dsensdt_vert_vert, 'dmoistdt':dmoistdt_vert,'dmoistdt_vert': dmoistdt_vert_vert, 'dpotdt':dpotdt_vert, 'dpotdt_vert':dpotdt_vert_vert}

logging.debug("loaded tendency dictionary")

save(destination+"tendency_Wm2_dic.hkl"      ,tendency_Wm2_dic)

logging.debug("successfully saved dictionaries in files")


################################################################ 




