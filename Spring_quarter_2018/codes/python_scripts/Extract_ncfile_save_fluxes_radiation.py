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
log_directory='/project2/tas1/pragallva/Spring_quarter_2018/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
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
num=str(dirc[3]) ## represents an annual year of data
source='/project2/tas1/pragallva/Spring_quarter_2018/exp_data/'+dirc[1]+'_'+dirc[2]+'/'
one_year=source+dirc[1]+'_'+dirc[2]+num+'.nc'

destination='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+num+'/'
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
EKE       =v_var['EKE'][...]  ### Specific humidity
logging.debug('saved EKE')

#virtual_T = (1.0+q*(Rv/Rd-1))*temp
#logging.debug('calculated virtual T')
# In[13]:


u_sq= u_comp**2
v_sq= v_comp**2
KE  = (u_sq+v_sq)/2.0

logging.debug('calculated KE')

CpT       =  Cp*temp
Lq        =  L*q
gZ        =  g*Z
# gZ_calc =  g*Z_calc
MSE       = (CpT+Lq+gZ+KE)
MSE_flux  = v_comp * (MSE)

logging.debug('calculated MSE fluxes')

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)


raw_data_dic      ={'Cp':Cp,"temp":R(temp),"u":R(u_comp),'v':R(v_comp),'L':L,'q':R(q),'g':g,'Z':R(Z) }
raw_MSE_dic       ={'MSE_flux':R(MSE_flux) }

logging.debug('saved raw data')

save(destination+"T_uv_Z.hkl",raw_data_dic)
save(destination+"all_MSE_flux.hkl",raw_MSE_dic)


def vertical_integral(arg,pressure):
    integral=0
    for i in range(0,len(sigma_full)):
        integral=integral+arg[:,i,:,:]*(-pressure[:,i+1,:,:]+pressure[:,i,:,:])/g
    return integral

vert_MSE_flux       = vertical_integral(MSE_flux  ,pres_half)
logging.debug('vertical integartion : MSE fluxes')
vert_sensible_flux  = vertical_integral(CpT*v_comp,pres_half)
logging.debug('vertical integration : sensible fluxes')
vert_latent_flux    = vertical_integral(Lq*v_comp ,pres_half)
logging.debug('vertical integration : latent fluxes')
vert_potential_flux = vertical_integral(gZ*v_comp ,pres_half)
logging.debug('vertical integration : nc potential fluxes')
#vert_potential_flux_calc = vertical_integral(gZ_calc*v_comp ,p_half)
#logging.debug('vertical integration : calc potential fluxes')
vert_KE_flux        = vertical_integral(KE*v_comp ,pres_half)
logging.debug('vertical integration : kinetic energy fluxes')
vert_enthalpy       = vertical_integral((CpT+Lq) ,pres_half)
logging.debug('vertical integration : Enthalpy ')


# RADIATION DATA

# In[ ]:

SW_sfc=v_var['flux_sw'][...] # Net SW surface flux
LW_sfc_dn=v_var['flux_lw'][...] # LW surface flux down
olr   =v_var['olr'][...]     # Outgoing LW radiation
SW_toa=v_var['toa_sw'][...]  # Net TOA SW flux

shflx=v_var['flux_t'][...]   ### Surface sensible heat flux
lhflx=v_var['flux_lhe'][...] #### Latent heat of fusion 
logging.debug('saved radiation')

t_surf    =v_var['t_surf'][...]  ### Specific humidity
logging.debug('saved t_surf')

flux_oceanq =v_var['flux_oceanq'][...]  ### Specific humidity
logging.debug('saved flux_oceanq')

LW_sfc_up=sigma*t_surf**4 # LW surface flux down


TOA= +SW_toa - olr                   ## downwards
SFC= shflx + lhflx + LW_sfc_up - LW_sfc_dn - SW_sfc   ## upwards
Net_rad=SFC+TOA

SWABS= +SW_toa - SW_sfc
SHF  = shflx + lhflx + LW_sfc_up - LW_sfc_dn
Net_rad2= SWABS + SHF - olr


logging.debug("radiation-- I think its right now !")


# $\frac{dm}{dt}$

# In[ ]:

#(1440, 64, 128)
dh_by_dt=np.copy(vert_enthalpy)
dt=6*60*60
for la in range(dh_by_dt.shape[1]):
    for lo in range(dh_by_dt.shape[-1]):
        dh_by_dt[:,la,lo]=np.gradient( vert_enthalpy[:,la,lo],dt)

logging.debug("calculated dh_by_dt")

# SAVING AS A DICTIONARY

plev = np.array([0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0])

coord_dic       ={"lat":lat,"lon":lon,"time":times,"p_sfc":p_sfc,"no_of_plevels":len(sigma_full),"plev": plev} 

logging.debug("loaded coordinates dictionary")

fluxes_dic      ={"MSE_flux":vert_MSE_flux,"sensible_flux":vert_sensible_flux,'latent_flux':vert_latent_flux,'potential_flux':vert_potential_flux,'KE_flux':vert_KE_flux,'dh_by_dt':dh_by_dt}

logging.debug("loaded MSE_flux dictionary")

#rad_dic         ={"SW_sfc_d":SW_sfc,"LW_sfc_d":LW_sfc_dn-LW_sfc_up,'olr':olr,'SW_toa_d':SW_toa,'shflx_u':shflx,'lhflx_u':lhflx,'TOA_d':TOA,'SFC_u':SFC,'Net_rad':Net_rad, 'SWABS':SWABS, 'SHF':SHF, 'flux_oceanq':flux_oceanq}

logging.debug("loaded rad_dic dictionary")

# In[141]:
save(destination+"coord_dic.hkl"       ,coord_dic)
np.save(destination+"time.npy"         ,time)
np.save(destination+"EKE.npy"          ,EKE)
save(destination+"fluxes_dic.hkl"      ,fluxes_dic)
#save(destination+"rad_dic.hkl"         ,rad_dic)

logging.debug("successfully saved dictionaries in files")


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

def MSE_total(m):
    v=v_comp
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    flux      = (V)*weighted(M)                # (year, month,days, hour, plev, lat, lon)
    monthly_f = flux.mean(axis=2).mean(axis=2) # (year, month, plev, lat, lon)
    zonal_f   = monthly_f.mean(axis=-1)        # (year, month, plev, lat)
    vert_flux = zonal_f.sum(axis=2)            # (year, month, lat)
    return vert_flux.mean(axis=0)              # (month, lat)--> yearly average

def mean_meridional(m):
    v=v_comp
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = (M).mean(axis=2).mean(axis=2)  # (year, month, days, hour, plev, lat, lon)
    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    zonal_m   = monthly_m.mean(axis=-1)        # (year, month, plev, lat)
    zonal_v   = monthly_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=2)  # (year, month, lat)
    return vert_flux.mean(axis=0)              # (month, lat)--> yearly average
    
def stationary_eddies(m):
    v=v_comp
    M         = reshape(m)                              # (year, month, days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = (M).mean(axis=2).mean(axis=2)           # (year, month,days, hour, plev, lat, lon)
    monthly_v = (V).mean(axis=2).mean(axis=2)                        # (year, month, plev, lat, lon)
    m_star    = monthly_m-monthly_m.mean(axis=-1)[...,None]          # (year, month, plev, lat, lon)
    v_star    = monthly_v-monthly_v.mean(axis=-1)[...,None]
    flux_weighted  = weighted((m_star*v_star)[:,:,None,None,:,:,:])    # (year, month, plev, lat)
    vert_flux = flux_weighted.mean(axis=2).mean(axis=2).mean(axis=-1).sum(axis=2)  # (year, month, lat)
    return vert_flux.mean(axis=0)              # (month, lat)--> yearly average
    
def transient_eddies(m):
    v=v_comp
    M         = reshape(m)                            # (year, month,days, hour, plev, lat, lon)
    V         = reshape(v)  
    monthly_m = M.mean(axis=2).mean(axis=2)           # average over days 
    monthly_v = V.mean(axis=2).mean(axis=2)           # (month, plev, lat, lon)
    m_prime   = M-monthly_m[:,:,None,None,:,:,:]      # (month, days, plev, lat, lon)
    v_prime   = V-monthly_v[:,:,None,None,:,:,:]
    flux      = weighted(m_prime*v_prime).mean(axis=2).mean(axis=2).mean(axis=-1) # (year,month,plev,lat) 
    vert_flux = (flux).sum(axis=2)                     # (month, year, lat) 
    return vert_flux.mean(axis=0)                      # (month, lat)--> yearly average

def eke_flux():
    v         =v_comp
    u         =u_comp
    U         = reshape(every_day_readjust(u))
    V         = reshape(every_day_readjust(v)) 
    KE        = (U**2+V**2)/2.0
    flux      = (KE*V).mean(axis=2).mean(axis=-1)               # (month,  plev, lat)   
    vert_flux = vertical_integral(flux)                         # (month,  lat) 
    return vert_flux

#############################
########   dm_by_dt    ######
#############################

def tendency_six_hours(h):
    dt=6*60*60
    dh_by_dt=np.copy(h)    
    for lo in range(len(lon)):
        for la in range(len(lat)):
            for lev in range(len(sigma_full)):
                dh_by_dt[:,lev,la,lo]=np.gradient( h[:,lev,la,lo],dt)                
    return dh_by_dt

logging.debug("calculated dh_by_dt dictionary")

moist_enthalpy = CpT+Lq+KE   ## ((1440, 29, 64, 128))
dh_by_dt=tendency_six_hours(moist_enthalpy)

dhdt=weighted(reshape(dh_by_dt)).mean(axis=2).mean(axis=2).mean(axis=-1) # reshape 
dhdt_vert=dhdt.sum(axis=2).mean(axis=0)              # (month, lat)--> yearly average and vertical average
# dhdt=reshape(every_day_readjust(dh_by_dt))

mse_flux=MSE_total(CpT+gZ+Lq+KE)

# In[48]:

MM_flux     =mean_meridional(CpT+gZ+Lq)
SE_flux     =stationary_eddies(CpT+gZ+Lq)
TE_flux     =transient_eddies(CpT+gZ+Lq)

KE_flux=MSE_total(KE)

TE_sensible=transient_eddies(CpT)
TE_moist   =transient_eddies(Lq)
TE_pot     =transient_eddies(gZ)

SE_sensible=stationary_eddies(CpT)
SE_moist   =stationary_eddies(Lq)
SE_pot     =stationary_eddies(gZ)

MM_sensible=mean_meridional(CpT)
MM_moist   =mean_meridional(Lq)
MM_pot     =mean_meridional(gZ)

# saving time as  raw format


zonal_decomposed_fluxes_dic      ={"MM_flux":MM_flux,"SE_flux":SE_flux,'TE_flux':TE_flux,'KE_flux':KE_flux,'TE_sensible':TE_sensible,'TE_moist':TE_moist,"SE_sensible":SE_sensible,"SE_moist":SE_moist,"SE_pot":SE_pot,"MM_sensible":MM_sensible,"MM_moist":MM_moist,"MM_pot":MM_pot,"dhdt":dhdt_vert}

save(destination+"zonal_decomposed_fluxes_dic.hkl" ,zonal_decomposed_fluxes_dic)
logging.debug("loaded zonal decomposed fluxes dictionary")

def reshape_rad(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(lat),len(lon)))
    return y1.mean(axis=-1).mean(axis=2).mean(axis=2).mean(axis=0)


zonal_radiation_dic         ={"SW_sfc_d":reshape_rad(SW_sfc),"LW_sfc_d":reshape_rad(LW_sfc_dn-LW_sfc_up),'olr':reshape_rad(olr),'SW_toa_d':reshape_rad(SW_toa),'shflx_u':reshape_rad(shflx),'lhflx_u':reshape_rad(lhflx),'TOA_d':reshape_rad(TOA),'SFC_u':reshape_rad(SFC),'Net_rad':reshape_rad(Net_rad), 'SWABS':reshape_rad(SWABS), 'SHF':reshape_rad(SHF), 'flux_oceanq':reshape_rad(flux_oceanq)}

save(destination+"zonal_radiation_dic.hkl",zonal_radiation_dic )
logging.debug("loaded zonal radiation dictionary")


# In[ ]:



# ####################
# #### soomthening ###
# ####################
# def smooth(y, box_pts):
#     box_pts=5
#     box = np.ones(box_pts)/box_pts
#     y_smooth=np.copy(y)
#     for m in range(mon):
#         y_smooth[m,:] = np.convolve(y[m,:], box, mode='same')
#     return y_smooth.transpose()

# print "successfully completed"

