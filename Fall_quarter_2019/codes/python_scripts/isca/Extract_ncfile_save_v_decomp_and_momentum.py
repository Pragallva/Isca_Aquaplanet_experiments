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

exp_dir=str(dirc[1])
num    =str(dirc[2])
no_of_days  =dirc[3]

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

log_directory = quarter+'/codes/shell_scripts/log/'
make_sure_path_exists( path   = log_directory )


for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)


logging.basicConfig( filename = log_directory+'momentum'+exp_dir+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

source   = quarter+'/exp_data/isca_repeat/'+exp_dir+'/'
one_year = source+exp_dir+num+'.nc'

destination = quarter+'/post_process_data/isca_repeat/avged_over'+str(no_of_days)+'days/'+exp_dir+num+'/'
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


if no_of_days == None:
   no_of_days = 30
else:
   days = int(no_of_days)

no_of_years = len(time)/(12*30*4)
no_of_months =12
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
omega     =v_var['omega'][:,::-1,:,:]
logging.debug('saved omega')
Z         =v_var['height'][:,::-1,:,:]
logging.debug('saved Z')
q         =v_var['sphum'][:,::-1,:,:]  ### Specific humidity
logging.debug('saved q')

u_sq= u_comp**2
v_sq= v_comp**2
KE  = (u_sq+v_sq)/2.0

logging.debug('calculated KE')

CpT       =  Cp*temp
Lq        =  L*q
gZ        =  G*Z
MSE       = (CpT+Lq+gZ+KE)

logging.debug('calculated MSE fluxes')

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)


logging.debug("loaded coordinates dictionary")


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
weights = (ps[:,:,:,:,:-1,:,:]-ps[:,:,:,:,1:,:,:])/G

def weighted(arg):
    w = arg*weights
    return w
    
############## decompose v into v_MM, v_TE and v_SE ###################

def zonal_mean(X):
    return np.nanmean(X, axis=-1)[..., None]

def time_mean(X):
    daily_mean   = np.nanmean(X, axis=3)[:,:,:,None,:,:,:]   ## Avged across all hours
    monthly_mean = np.nanmean(daily_mean, axis=2)[:,:,None,:,:,:,:] ## Avged across all months
    return monthly_mean

def time_eddy(X):
    return X-time_mean(X)

def zonal_eddy(X):
    return X-time_mean(X)

def zonal_eddy_time_mean_flux(g=u_comp, h=v_comp): #### This is similar to stationary eddies as defined in our paper
    A                       = reshape(g)
    B                       = reshape(h)
    zonal_eddy_A            = zonal_eddy(A)
    zonal_eddy_B            = zonal_eddy(B)
    time_mean_zonal_eddy_A  = time_mean(zonal_eddy_A)
    time_mean_zonal_eddy_B  = time_mean(zonal_eddy_B)
    flux                    = time_mean_zonal_eddy_A*time_mean_zonal_eddy_A
    return zonal_mean(flux)

def zonal_mean_time_mean_flux(g=u_comp, h=v_comp):
    A                       = reshape(g)     # (year, month,days, hour, plev, lat, lon)
    B                       = reshape(h)
    zonal_A                 = zonal_mean(A)
    zonal_B                 = zonal_mean(B)
    time_zonal_mean_A       = time_mean(zonal_A)
    time_zonal_mean_B       = time_mean(zonal_B)
    flux                    = time_zonal_mean_A*time_zonal_mean_B
    return flux

def time_eddy_flux(g=u_comp, h=v_comp): # This is similar to transient eddies as defined in our paper
    A                       = reshape(g)
    B                       = reshape(h)
    time_eddy_A             = time_eddy(A)
    time_eddy_B             = time_eddy(B)
    flux                    = time_eddy_A*time_eddy_B
    return zonal_mean(time_mean(flux))

def dteddy_dp(g=u_comp, h=v_comp, ax=-3):
    dp                      = zonal_mean(time_mean(weights*G))
    A                       = reshape(g)
    B                       = reshape(h)
    time_eddy_A             = time_eddy(A)
    time_eddy_B             = time_eddy(B)
    flux                    = time_eddy_A*time_eddy_B
    Flux                    = np.append(flux, flux[...,-1,None,:,:],axis=-3)
    dflux                   = zonal_mean(time_mean(Flux[...,:-1,:,:]-Flux[...,1:,:,:]))
    dudp                    = dflux/dp
    return dudp

def dzeddy_tmean_dp(g=u_comp, h=v_comp, ax=-3):
    dp                      = zonal_mean(time_mean(weights*G))
    A                       = reshape(g)
    B                       = reshape(h)
    zonal_eddy_A            = zonal_eddy(A)
    zonal_eddy_B            = zonal_eddy(B)
    time_mean_zonal_eddy_A  = time_mean(zonal_eddy_A)
    time_mean_zonal_eddy_B  = time_mean(zonal_eddy_B)
    flux                    = time_mean_zonal_eddy_A*time_mean_zonal_eddy_A
    Flux                    = np.append(flux, flux[...,-1,None,:,:],axis=-3)
    dflux                   = zonal_mean(Flux[...,:-1,:,:]-Flux[...,1:,:,:])
    dudp                    = dflux/dp
    return dudp


def vort(g=u_comp,ax=2):
   radius       =  a
   cos_phi      =  np.cos(np.deg2rad(lat))[None, None, :, None]
   vorticity    = -np.gradient(g*cos_phi, np.deg2rad(lat), axis = ax)/(radius*cos_phi)
   return vorticity

def ugrad_dp(g=u_comp, ax=1):
    dp                      = pres_half[:,:-1,...]-pres_half[:,1:,...]
    g                       = np.append(g, g[:,-1,None,...],axis=1) 
    du                      = g[:,:-1,...]-g[:,1:,...]
    dudp                    = du/dp
    #DuDp                   = reshape(dudp)
    return dudp

def u_tendency(g=u_comp):
    dt                      = 6*60*60
    dudt                    = np.gradient(g, dt, axis=0)
#    DuDt                    = reshape(dudt)
    return dudt

def div_EMF(X):
  radius       = a
  cos_phi      = (np.cos(np.deg2rad(lat)))[None,None,None,None,None,:,None]
  div = (1/(radius*cos_phi**2)) * (np.gradient (X*cos_phi**2, np.deg2rad(lat), axis=-2 ))
  return div

rot_rate         = 7.29*10**(-5)
sin_phi          = (np.sin(np.deg2rad(lat)))
f                = (2*rot_rate*(sin_phi))
F                = f[None,None,None,None,None,:,None]

def ZTM(Y):
   return np.squeeze(zonal_mean(time_mean(Y)))

zeta_mean        = (zonal_mean(time_mean(reshape(vort())) ))  
logging.debug(str(num)+'calculated zeta')

du_dt            =  (zonal_mean(time_mean(reshape(u_tendency())) ))
logging.debug(str(num)+'calculated dudt, shape = '+ str(np.shape(du_dt)))

w_dudp_zmean_tmean = zonal_mean_time_mean_flux( g=omega, h=ugrad_dp() )
logging.debug(str(num)+'calculated w_dudp, shape = '+ str(np.shape(w_dudp_zmean_tmean)))

zeta_v_zmean_tmean =  zonal_mean_time_mean_flux( g=vort(), h=v_comp )
logging.debug(str(num)+'calculated zeta_v, shape = '+ str(np.shape(zeta_v_zmean_tmean)))

######

div_EMF_teddy   = div_EMF(time_eddy_flux(g=u_comp, h=v_comp))
logging.debug(str(num)+'calculated div TE EMF, shape = '+ str(np.shape(div_EMF_teddy)))

du_w_teddy_dp   =  dteddy_dp(g=u_comp, h=omega, ax=-3)
logging.debug(str(num)+"calculated TE d_u'w'_dp, shape = "+ str(np.shape(du_w_teddy_dp)))

#####


div_EMF_zeddy_tmean   = div_EMF(zonal_eddy_time_mean_flux(g=u_comp, h=v_comp))
logging.debug(str(num)+'calculated div SE EMF, shape = '+ str(np.shape(div_EMF_zeddy_tmean)))

du_w_zeddy_tmean_dp   =  dzeddy_tmean_dp(g=u_comp, h=omega, ax=-3)
logging.debug(str(num)+"calculated SE d_u'w'_dp, shape = "+ str(np.shape(du_w_zeddy_tmean_dp)))

#####

def ZT(Y):
   return (zonal_mean(time_mean(Y)))

v_dudt  = du_dt/F
v_zetaV = -zeta_v_zmean_tmean/F
v_wdudp = w_dudp_zmean_tmean/F
v_MM    = v_dudt + v_zetaV + v_wdudp

v_duw_dp_teddy  = du_w_teddy_dp/F
v_div_EMF_teddy = div_EMF_teddy/F
v_TE            = v_duw_dp_teddy + v_div_EMF_teddy

v_duw_dp_zeddy_tmean  = du_w_zeddy_tmean_dp/F
v_div_EMF_zeddy_tmean = div_EMF_zeddy_tmean/F
v_SE                  = v_duw_dp_zeddy_tmean + v_div_EMF_zeddy_tmean

v_net                 = ZT(reshape(v_comp))-(ZT(v_TE)+ZT(v_SE)+ZT(v_MM))

logging.debug(str(num)+'calculated all vs')

####

mom_budget = {'du_dt':ZTM(du_dt),'zeta_v_zmean_tmean':ZTM(zeta_v_zmean_tmean), \
              'w_dudp_zmean_tmean':ZTM(w_dudp_zmean_tmean), 'du_w_teddy_dp':ZTM(du_w_teddy_dp), \
              'du_w_zeddy_tmean_dp':ZTM(du_w_zeddy_tmean_dp), 'div_EMF_teddy':ZTM(div_EMF_teddy), \
              'div_EMF_zeddy_tmean':ZTM(div_EMF_zeddy_tmean)}

logging.debug(str(num)+'Saved momentum budget dictionary')
#####

v_dudt1  = du_dt/(F+zeta_mean)
v_wdudp1 = w_dudp_zmean_tmean/(F+zeta_mean)
v_MM1    = v_dudt1 + v_wdudp1

v_duw_dp_teddy1  = du_w_teddy_dp/(F+zeta_mean)
v_div_EMF_teddy1 = div_EMF_teddy/(F+zeta_mean)
v_TE1            = v_duw_dp_teddy1 + v_div_EMF_teddy1

v_duw_dp_zeddy_tmean1  = du_w_zeddy_tmean_dp/(F+zeta_mean)
v_div_EMF_zeddy_tmean1 = div_EMF_zeddy_tmean/(F+zeta_mean)
v_SE1                  = v_duw_dp_zeddy_tmean1 + v_div_EMF_zeddy_tmean1

v_net1                 = ZT(reshape(v_comp))-(ZT(v_TE1)+ZT(v_SE1)+ZT(v_MM1))

logging.debug(str(num)+'calculated all v1 s')

####

vs = {'v_tot':ZTM(reshape(v_comp)),'v_TE':ZTM(v_TE), 'v_SE':ZTM(v_SE), 'v_MM':ZTM(v_MM), \
      'v_TE1':ZTM(v_TE1), 'v_SE1':ZTM(v_SE1), 'v_MM1':ZTM(v_MM1),'v_net':ZTM(v_net),'v_net1':ZTM(v_net1),\
      'f':ZTM(F),'zeta':ZTM(zeta_mean)}

logging.debug(str(num)+'Saved v dictionary')

########### v decomposed mean meridional MSE flux   ############

def mean_meridional(m,vdecomp):
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = vdecomp
    monthly_m = (M).mean(axis=2).mean(axis=2)  # (year, month, days, hour, plev, lat, lon)
    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    zonal_m   = monthly_m.mean(axis=-1)        # (year, month, plev, lat)
    zonal_v   = monthly_v.mean(axis=-1)
    vert_flux = (zonal_m*zonal_v).sum(axis=2)  # (year, month, lat)
    return vert_flux.mean(axis=0),  (zonal_m*zonal_v).mean(axis=0)


v_MM_flux, v_MM_flux_vert = mean_meridional(CpT+gZ+Lq, v_MM)
v_SE_flux, v_SE_flux_vert = mean_meridional(CpT+gZ+Lq, v_SE)
v_TE_flux, v_TE_flux_vert = mean_meridional(CpT+gZ+Lq, v_TE)
v_net_flux, v_net_flux_vert = mean_meridional(CpT+gZ+Lq, v_net)
v_MM1_flux, v_MM1_flux_vert = mean_meridional(CpT+gZ+Lq, v_MM1)
v_SE1_flux, v_SE1_flux_vert = mean_meridional(CpT+gZ+Lq, v_SE1)
v_TE1_flux, v_TE1_flux_vert = mean_meridional(CpT+gZ+Lq, v_TE1)
v_net1_flux, v_net1_flux_vert = mean_meridional(CpT+gZ+Lq, v_net1)
logging.debug(str(num)+'Calulated MMC flux')

v_MM_sens_flux, v_MM_sens_flux_vert = mean_meridional(CpT, v_MM)
v_SE_sens_flux, v_SE_sens_flux_vert = mean_meridional(CpT, v_SE)
v_TE_sens_flux, v_TE_sens_flux_vert = mean_meridional(CpT, v_TE)
v_net_sens_flux, v_net_sens_flux_vert = mean_meridional(CpT, v_net)
v_MM1_sens_flux, v_MM1_sens_flux_vert = mean_meridional(CpT, v_MM1)
v_SE1_sens_flux, v_SE1_sens_flux_vert = mean_meridional(CpT, v_SE1)
v_TE1_sens_flux, v_TE1_sens_flux_vert = mean_meridional(CpT, v_TE1)
v_net1_sens_flux, v_net1_sens_flux_vert = mean_meridional(CpT, v_net1)
logging.debug(str(num)+'Calulated MMC sens flux')

v_MM_pot_flux, v_MM_pot_flux_vert = mean_meridional(gZ, v_MM)
v_SE_pot_flux, v_SE_pot_flux_vert = mean_meridional(gZ, v_SE)
v_TE_pot_flux, v_TE_pot_flux_vert = mean_meridional(gZ, v_TE)
v_net_pot_flux, v_net_pot_flux_vert = mean_meridional(gZ, v_net)
v_MM1_pot_flux, v_MM1_pot_flux_vert = mean_meridional(gZ, v_MM1)
v_SE1_pot_flux, v_SE1_pot_flux_vert = mean_meridional(gZ, v_SE1)
v_TE1_pot_flux, v_TE1_pot_flux_vert = mean_meridional(gZ, v_TE1)
v_net1_pot_flux, v_net1_pot_flux_vert = mean_meridional(gZ, v_net1)
logging.debug(str(num)+'Calulated MMC potential flux')

v_MM_moist_flux, v_MM_moist_flux_vert = mean_meridional(Lq, v_MM)
v_SE_moist_flux, v_SE_moist_flux_vert = mean_meridional(Lq, v_SE)
v_TE_moist_flux, v_TE_moist_flux_vert = mean_meridional(Lq, v_TE)
v_net_moist_flux, v_net_moist_flux_vert = mean_meridional(Lq, v_net)
v_MM1_moist_flux, v_MM1_moist_flux_vert = mean_meridional(Lq, v_MM1)
v_SE1_moist_flux, v_SE1_moist_flux_vert = mean_meridional(Lq, v_SE1)
v_TE1_moist_flux, v_TE1_moist_flux_vert = mean_meridional(Lq, v_TE1)
v_net1_moist_flux, v_net1_moist_flux_vert = mean_meridional(Lq, v_net1)
logging.debug(str(num)+'Calulated MMC moist flux')

########### decomposing thermodynamic and dynamic parts  ##########


def delta_v_mmc(m,vdecomp):
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = vdecomp  
    monthly_m = (M).mean(axis=2).mean(axis=2)  # (year, month, days, hour, plev, lat, lon)
    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    zonal_m   = monthly_m.mean(axis=-1)        # (year, month, plev, lat)
    zonal_v   = monthly_v.mean(axis=-1)
    delta_v   = zonal_v-zonal_v.mean(axis=1)[:,None,...]
    vert_flux = (zonal_m*delta_v).sum(axis=2)  # (year, month, lat)
    return vert_flux.mean(axis=0), (zonal_m*delta_v).mean(axis=0)   
    
def delta_m_mmc(m,vdecomp):
    M         = reshape(m)                     # (year, month,days, hour, plev, lat, lon)
    V         = vdecomp  
    monthly_m = (M).mean(axis=2).mean(axis=2)  # (year, month, days, hour, plev, lat, lon)
    monthly_v = weighted(V).mean(axis=2).mean(axis=2)  # (year, month, plev, lat, lon)
    zonal_m   = monthly_m.mean(axis=-1)        # (year, month, plev, lat)
    delta_m   = zonal_m-zonal_m.mean(axis=1)[:,None,...]    
    zonal_v   = monthly_v.mean(axis=-1)
    vert_flux = (delta_m*zonal_v).sum(axis=2)  # (year, month, lat)
    return vert_flux.mean(axis=0), (delta_m*zonal_v).mean(axis=0)


v_MM_del_m_mmc, v_MM_del_m_mmc_vert    = delta_m_mmc(CpT+gZ+Lq, v_MM)
v_SE_del_m_mmc, v_SE_del_m_mmc_vert    = delta_m_mmc(CpT+gZ+Lq, v_SE)
v_TE_del_m_mmc, v_TE_del_m_mmc_vert    = delta_m_mmc(CpT+gZ+Lq, v_TE)
v_net_del_m_mmc, v_net_del_m_mmc_vert  = delta_m_mmc(CpT+gZ+Lq, v_net)
v_MM1_del_m_mmc, v_MM1_del_m_mmc_vert  = delta_m_mmc(CpT+gZ+Lq, v_MM1)
v_SE1_del_m_mmc, v_SE1_del_m_mmc_vert  = delta_m_mmc(CpT+gZ+Lq, v_SE1)
v_TE1_del_m_mmc, v_TE1_del_m_mmc_vert  = delta_m_mmc(CpT+gZ+Lq, v_TE1)
v_net1_del_m_mmc, v_net1_del_m_mmc_vert  = delta_m_mmc(CpT+gZ+Lq, v_net1)
logging.debug(str(num)+'Calulated MMC del m mmc flux')

v_MM_del_v_mmc, v_MM_del_v_mmc_vert       = delta_v_mmc(CpT+gZ+Lq, v_MM)
v_SE_del_v_mmc, v_SE_del_v_mmc_vert       = delta_v_mmc(CpT+gZ+Lq, v_SE)
v_TE_del_v_mmc, v_TE_del_v_mmc_vert       = delta_v_mmc(CpT+gZ+Lq, v_TE)
v_net_del_v_mmc, v_net_del_v_mmc_vert     = delta_v_mmc(CpT+gZ+Lq, v_net)
v_MM1_del_v_mmc, v_MM1_del_v_mmc_vert   = delta_v_mmc(CpT+gZ+Lq, v_MM1)
v_SE1_del_v_mmc, v_SE1_del_v_mmc_vert   = delta_v_mmc(CpT+gZ+Lq, v_SE1)
v_TE1_del_v_mmc, v_TE1_del_v_mmc_vert   = delta_v_mmc(CpT+gZ+Lq, v_TE1)
v_net1_del_v_mmc, v_net1_del_v_mmc_vert   = delta_v_mmc(CpT+gZ+Lq, v_net1)
logging.debug(str(num)+'Calulated MMC del v mmc flux')

v_MM_MMC = {"v_MM_flux":v_MM_flux, "v_MM_flux_vert":v_MM_flux_vert, \
            "v_MM_sens_flux": v_MM_sens_flux,"v_MM_sens_flux_vert":v_MM_sens_flux_vert,\
            "v_MM_pot_flux":v_MM_pot_flux,"v_MM_pot_flux_vert":v_MM_pot_flux_vert,\
            "v_MM_moist_flux":v_MM_moist_flux,"v_MM_moist_flux_vert":v_MM_moist_flux_vert,\
            "v_MM_del_m_mmc":v_MM_del_m_mmc,"v_MM_del_m_mmc_vert":v_MM_del_m_mmc_vert,\
            "v_MM_del_v_mmc":v_MM_del_v_mmc,"v_MM_del_v_mmc_vert":v_MM_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_MM_MMC dictionary')

v_MM1_MMC = {"v_MM1_flux":v_MM1_flux, "v_MM1_flux_vert":v_MM1_flux_vert, \
            "v_MM1_sens_flux": v_MM1_sens_flux,"v_MM1_sens_flux_vert":v_MM1_sens_flux_vert,\
            "v_MM1_pot_flux":v_MM1_pot_flux,"v_MM1_pot_flux_vert":v_MM1_pot_flux_vert,\
            "v_MM1_moist_flux":v_MM1_moist_flux,"v_MM1_moist_flux_vert":v_MM1_moist_flux_vert,\
            "v_MM1_del_m_mmc":v_MM1_del_m_mmc,"v_MM1_del_m_mmc_vert":v_MM1_del_m_mmc_vert,\
            "v_MM1_del_v_mmc":v_MM1_del_v_mmc,"v_MM1_del_v_mmc_vert":v_MM1_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_MM1_MMC dictionary')

v_TE_MMC = {"v_TE_flux":v_TE_flux, "v_TE_flux_vert":v_TE_flux_vert,\
            "v_TE_sens_flux": v_TE_sens_flux,"v_TE_sens_flux_vert":v_TE_sens_flux_vert,\
            "v_TE_pot_flux":v_TE_pot_flux,"v_TE_pot_flux_vert":v_TE_pot_flux_vert,\
            "v_TE_moist_flux":v_TE_moist_flux,"v_TE_moist_flux_vert":v_TE_moist_flux_vert,\
            "v_TE_del_m_mmc":v_TE_del_m_mmc,"v_TE_del_m_mmc_vert":v_TE_del_m_mmc_vert,\
            "v_TE_del_v_mmc":v_TE_del_v_mmc,"v_TE_del_v_mmc_vert":v_TE_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_TE_MMC dictionary')

v_TE1_MMC = {"v_TE1_flux":v_TE1_flux, "v_TE1_flux_vert":v_TE1_flux_vert,\
            "v_TE1_sens_flux": v_TE1_sens_flux,"v_TE1_sens_flux_vert":v_TE1_sens_flux_vert,\
            "v_TE1_pot_flux":v_TE1_pot_flux,"v_TE1_pot_flux_vert":v_TE1_pot_flux_vert,\
            "v_TE1_moist_flux":v_TE1_moist_flux,"v_TE1_moist_flux_vert":v_TE1_moist_flux_vert,\
            "v_TE1_del_m_mmc":v_TE1_del_m_mmc,"v_TE1_del_m_mmc_vert":v_TE1_del_m_mmc_vert,\
            "v_TE1_del_v_mmc":v_TE1_del_v_mmc,"v_TE1_del_v_mmc_vert":v_TE1_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_TE1_MMC dictionary')

v_SE_MMC = {"v_SE_flux":v_SE_flux, "v_SE_flux_vert":v_SE_flux_vert,\
            "v_SE_sens_flux": v_SE_sens_flux,"v_SE_sens_flux_vert":v_SE_sens_flux_vert,\
            "v_SE_pot_flux":v_SE_pot_flux,"v_SE_pot_flux_vert":v_SE_pot_flux_vert,\
            "v_SE_moist_flux":v_SE_moist_flux,"v_SE_moist_flux_vert":v_SE_moist_flux_vert,\
            "v_SE_del_m_mmc":v_SE_del_m_mmc,"v_SE_del_m_mmc_vert":v_SE_del_m_mmc_vert,\
            "v_SE_del_v_mmc":v_SE_del_v_mmc,"v_SE_del_v_mmc_vert":v_SE_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_SE_MMC dictionary')


v_SE1_MMC = {"v_SE1_flux":v_SE1_flux, "v_SE1_flux_vert":v_SE1_flux_vert,\
            "v_SE1_sens_flux": v_SE1_sens_flux,"v_SE1_sens_flux_vert":v_SE1_sens_flux_vert,\
            "v_SE1_pot_flux":v_SE1_pot_flux,"v_SE1_pot_flux_vert":v_SE1_pot_flux_vert,\
            "v_SE1_moist_flux":v_SE1_moist_flux,"v_SE1_moist_flux_vert":v_SE1_moist_flux_vert,\
            "v_SE1_del_m_mmc":v_SE1_del_m_mmc,"v_SE1_del_m_mmc_vert":v_SE1_del_m_mmc_vert,\
            "v_SE1_del_v_mmc":v_SE1_del_v_mmc,"v_SE1_del_v_mmc_vert":v_SE1_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_SE1_MMC dictionary')


v_net_MMC = {"v_net_flux":v_net_flux, "v_net_flux_vert":v_net_flux_vert,\
            "v_net_sens_flux": v_net_sens_flux,"v_net_sens_flux_vert":v_net_sens_flux_vert,\
            "v_net_pot_flux":v_net_pot_flux,"v_net_pot_flux_vert":v_net_pot_flux_vert,\
            "v_net_moist_flux":v_net_moist_flux,"v_net_moist_flux_vert":v_net_moist_flux_vert,\
            "v_net_del_m_mmc":v_net_del_m_mmc,"v_net_del_m_mmc_vert":v_net_del_m_mmc_vert,\
            "v_net_del_v_mmc":v_net_del_v_mmc,"v_net_del_v_mmc_vert":v_net_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_net_MMC dictionary')

v_net1_MMC = {"v_net1_flux":v_net1_flux, "v_net1_flux_vert":v_net1_flux_vert,\
              "v_net1_sens_flux": v_net1_sens_flux,"v_net1_sens_flux_vert":v_net1_sens_flux_vert,\
              "v_net1_pot_flux":v_net1_pot_flux,"v_net1_pot_flux_vert":v_net1_pot_flux_vert,\
              "v_net1_moist_flux":v_net1_moist_flux,"v_net1_moist_flux_vert":v_net1_moist_flux_vert,\
              "v_net1_del_m_mmc":v_net1_del_m_mmc,"v_net1_del_m_mmc_vert":v_net1_del_m_mmc_vert,\
              "v_net1_del_v_mmc":v_net1_del_v_mmc,"v_net1_del_v_mmc_vert":v_net1_del_v_mmc_vert}

logging.debug(str(num)+'Saved v_net1_MMC dictionary')

v_MMC = {"v_MM_MMC":v_MM_MMC, "v_TE_MMC":v_TE_MMC, "v_SE_MMC":v_SE_MMC,\
         "v_MM1_MMC":v_MM1_MMC, "v_TE1_MMC":v_TE1_MMC, "v_SE1_MMC":v_SE1_MMC,\
         "v_net_MMC":v_net_MMC, "v_net1_MMC":v_net1_MMC}

logging.debug(str(num)+'Saved v_MMC dictionary')


save(destination+"v_MMC.hkl",v_MMC)
logging.debug("loaded MMC decomposed dictionary")

save(destination+"vs.hkl",vs)
logging.debug("loaded vs decomposed dictionary")

save(destination+"mom_budget.hkl",vs)
logging.debug("loaded momentum budget dictionary")

end = ti.time()
logging.debug("Time taken --> "+ str(end-start))
logging.debug("Awesome! complete!")
logging.debug("-------------------------------")

