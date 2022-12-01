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

import time as ti
start = ti.time()

exp_dir    = str(dirc[1]) ## "HC0_la50m_oc50m"
num        = str(dirc[2]) ## represents an annual year of data
no_of_days = dirc[3]


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

logging.basicConfig( filename = log_directory+'eddy'+exp_dir+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

source   = quarter+'/exp_data/isca_repeat/data_in_pres_coord/'+exp_dir+'/'
one_year = source+exp_dir+num+'.nc'

destination = quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+num+'/'
make_sure_path_exists(destination);

v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('....Ssave_MSE_fluxes....'+str(num)+' averaged over '+str(no_of_days)+' days')

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
rot_rate = 7.29*10**(-5) ## radian/s


def get_data(variable='transient', f = one_year):
    ncfile = f
    v_var  = nc.Dataset(ncfile,'r')
    data   = v_var.variables
    each_model_name = data[variable][:]
    v_var.close()
    return each_model_name

ucomp    = get_data('ucomp')
logging.debug(str(num)+'imported u')

vcomp    = get_data('vcomp')
logging.debug(str(num)+'imported v')

omega    = get_data('omega')
logging.debug(str(num)+'imported w')

height   = get_data('height')
logging.debug(str(num)+'imported Z')

sphum    = get_data('sphum')
logging.debug(str(num)+'imported q')

temp     = get_data('temp')
logging.debug(str(num)+'imported temp')

lat      = get_data('lat')
lon      = get_data('lon')
pres     = get_data('level')
time     = get_data('time')

logging.debug(str(num)+'imported coordinates')

coord = {'lat':lat, 'lon':lon, 'pres':pres, 'time':time}

if no_of_days == None:
   no_of_days = 30
else:
   days = int(no_of_days)

no_of_months = 360/days
hours        = 4
no_of_years  = len(time)/(no_of_months*days*hours)

u_sq  = ucomp**2
v_sq  = vcomp**2
KE  = (u_sq+v_sq)/2.0
#uw  = ucomp*omega
#uv  = ucomp*vcomp

cos_phi =  np.cos(np.deg2rad(lat))[None, None, :, None]
vort    = -np.gradient(ucomp * cos_phi, np.deg2rad(lat[1]-lat[0]), axis = -2)/(a*cos_phi)

logging.debug(str(num)+'calculated covariances, uu, vv, uw, uv, KE')

#CpT       =  Cp*temp
#Lq        =  L*q
#gZ        =  g*Z
#MSE       =  CpT+Lq+gZ+KE
#MSE_flux  =  v_comp * (MSE)

logging.debug(str(num)+'calculated MSE fluxes')

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(pres),len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)

def reshape(y):
    y1 = y.reshape((no_of_years,no_of_months,days,hours,len(pres),len(lat),len(lon)))
    return y1

#raw_data_dic    =  {'omega':R(omega),'Cp':Cp,"temp":R(temp),"u":R(ucomp),'v':R(vcomp),'L':L,'q':R(q),'g':g,'Z':R(Z) }
#raw_MSE_dic     =  {'MSE_flux':R(MSE_flux) }

#logging.debug(str(num)+'saved raw data')

##save(destination+"T_UV_Z.hkl",raw_data_dic)
##save(destination+"MSE_flux.hkl",raw_MSE_dic)

def zonal_mean(X):
    return np.nanmean(X, axis=-1)[..., None]

def time_mean(X):
    print X.shape
    daily_mean   = np.nanmean(X, axis=3)[:,:,:,None,:,:,:]          ## Averaged across all hours
    monthly_mean = np.nanmean(daily_mean, axis=2)[:,:,None,:,:,:,:] ## Averaged across all months
    return monthly_mean

def time_eddy(X):
    return X-time_mean(X)

def zonal_eddy(X):
    return X-time_mean(X)

def zonal_mean_time_mean_flux(a=vort, b=vcomp):    
    A                       = reshape(a)                     # (year, month,days, hour, plev, lat, lon)
    B                       = reshape(b)
    zonal_A                 = zonal_mean(A)
    zonal_B   	            = zonal_mean(B)
    time_zonal_mean_A       = time_mean(zonal_A)            
    time_zonal_mean_B       = time_mean(zonal_B)
    flux                    = time_zonal_mean_A*time_zonal_mean_B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def zonal_mean_time_eddy_flux(a=vort, b=vcomp):
    A                       = reshape(a)                     # (year, month,days, hour, plev, lat, lon)
    B                       = reshape(b)
    zonal_A                 = zonal_mean(A)
    zonal_B                 = zonal_mean(B)
    time_eddy_zonal_mean_A  = time_eddy(zonal_A)
    time_eddy_zonal_mean_B  = time_eddy(zonal_B)                
    flux                    = time_eddy_zonal_mean_A*time_eddy_zonal_mean_B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def zonal_eddy_time_eddy_flux(a=vort, b=vcomp):
    A                       = reshape(a)
    B                       = reshape(b)
    zonal_eddy_A            = zonal_eddy(A)
    zonal_eddy_B            = zonal_eddy(B)
    time_eddy_zonal_eddy_A  = time_eddy(zonal_eddy_A)
    time_eddy_zonal_eddy_B  = time_eddy(zonal_eddy_B)
    flux                    = time_eddy_zonal_eddy_A*time_eddy_zonal_eddy_B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def zonal_eddy_time_mean_flux(a=vort, b=vcomp): #### This is similar to stationary eddies as defined in our paper
    A                       = reshape(a)
    B                       = reshape(b)
    zonal_eddy_A            = zonal_eddy(A)
    zonal_eddy_B            = zonal_eddy(B)
    time_mean_zonal_eddy_A  = time_mean(zonal_eddy_A)
    time_mean_zonal_eddy_B  = time_mean(zonal_eddy_B)
    flux                    = time_mean_zonal_eddy_A*time_mean_zonal_eddy_A
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def zonal_mean_flux(a=vort, b=vcomp):
    A                       = reshape(a)
    B                       = reshape(b)
    zonal_mean_A            = zonal_mean(A)
    zonal_mean_B            = zonal_mean(B)
    flux                    = zonal_mean_A*zonal_mean_B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def zonal_eddy_flux(a=vort, b=vcomp):
    A                       = reshape(a)
    B                       = reshape(b)
    zonal_eddy_A            = zonal_eddy(A)
    zonal_eddy_B            = zonal_eddy(B)
    flux                    = zonal_eddy_A*zonal_eddy_B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def zonal_eddy_flux2(a=vort, b=vcomp):  ### mathematically this should be same as zonal_eddy_flux()
    x1                      = zonal_eddy_time_mean_flux(a,b)
    x2                      = zonal_eddy_time_eddy_flux(a,b)
    return x+x2

def time_eddy_flux(a=vort, b=vcomp):    #### This is similar to transient eddies as defined in our paper
    A                       = reshape(a)
    B                       = reshape(b)
    time_eddy_A             = time_eddy(A)
    time_eddy_B             = time_eddy(B)
    flux                    = time_eddy_A*time_eddy_B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

def time_eddy_flux2(a=vort, b=vcomp):  ### mathematically this should be same as zonal_eddy_flux()
    x1                      = zonal_mean_time_eddy_flux(a,b)
    x2                      = zonal_eddy_time_eddy_flux(a,b)
    return x+x2

def total_flux(a=vort, b=vcomp):
    A                       = reshape(a)
    B                       = reshape(b)
    flux                    = A*B
    mean_flux               = zonal_mean(time_mean(flux))
    return np.squeeze(mean_flux)

######## Now let us calculate the fluxes necessary for momentum fluxes ##########

def u_tendency(a=ucomp):
    dt                      = 6*60*60
    dudt                    = np.gradient(a, dt, axis=0)
    DuDt                    = reshape(dudt)
    print DuDt.shape
    return DuDt


def ugrad_dp(a=ucomp, ax=1):
    p                       = pres
    dudp                    = np.gradient(a, p, axis=ax)
    #DuDp                    = reshape(dudp)
    return dudp

def vort(a=ucomp):
   cos_phi      =  np.cos(np.deg2rad(lat))[None, None, :, None]
   vorticity    = -np.gradient(a * cos_phi, np.deg2rad(lat[1]-lat[0]), axis = -2)/(a*cos_phi)
#  VORTICITY    = reshape(VORTICITY)
   return vorticity

du_dt              = np.squeeze(zonal_mean(time_mean( u_tendency() )))
logging.debug(str(num)+'calculated dudt')



w_dudp_zmean_tmean = zonal_mean_time_mean_flux( a=omega, b=ugrad_dp() ) 
logging.debug(str(num)+'calculated w_dudp1')

w_dudp_zmean_teddy = zonal_mean_time_eddy_flux( a=omega, b=ugrad_dp() ) ### Zonal mean transient eddy fluxes
logging.debug(str(num)+'calculated w_dudp2')

w_dudp_zmean       = zonal_mean_flux(a=omega, b=ugrad_dp())
logging.debug(str(num)+'calculated w_dudp3')

w_dudp_zmean2      = w_dudp_zmean_tmean + w_dudp_zmean_teddy
logging.debug(str(num)+'calculated w_dudp4')

w_dudp = {'w_dudp_zmean_tmean': w_dudp_zmean_tmean, 'w_dudp_zmean_teddy': w_dudp_zmean_teddy, \
          'w_dudp_zmean': w_dudp_zmean, 'w_dudp_zmean2': w_dudp_zmean2}

# w_dudp_zeddy_teddy = zonal_eddy_time_eddy_flux( a=omega, b=ugrad_dp() ) ### Zonal mean transient eddy fluxes
# logging.debug(str(num)+'calculated w_dudp5')



zeta_v_zmean_tmean =  zonal_mean_time_mean_flux( a=vort(), b=vcomp )
logging.debug(str(num)+'calculated zeta_v1')

zeta_v_zmean_teddy =  zonal_mean_time_eddy_flux( a=vort(), b=vcomp )
logging.debug(str(num)+'calculated zeta_v2')

zeta_v_zmean       =  zonal_mean_flux( a=vort(), b=vcomp )
logging.debug(str(num)+'calculated zeta_v3')

zeta_v_zmean2      =  zeta_v_zmean_tmean + zeta_v_zmean_teddy
logging.debug(str(num)+'calculated zeta_v4')

zeta_v = {'zeta_v_zmean_tmean': zeta_v_zmean_tmean, 'zeta_v_zmean_teddy': zeta_v_zmean_teddy, \
          'zeta_v_zmean': zeta_v_zmean, 'zeta_v_zmean2': zeta_v_zmean2}


sin_phi          = np.deg2rad(np.sin(lat))[None,None,:,None]
fv               = 2*rot_rate*sin_phi*reshape(vcomp)
fv_zmean_tmean   = np.squeeze(zonal_mean(time_mean(fv)))


EMF_zeddy        =  zonal_eddy_flux(a=ucomp, b=vcomp)
logging.debug(str(num)+'calculated EMF1')

EMF_zeddy_tmean  =  zonal_eddy_time_mean_flux(a=ucomp, b=vcomp)
logging.debug(str(num)+'calculated EMF2')

EMF_zeddy_teddy  =  zonal_eddy_time_eddy_flux(a=ucomp, b=vcomp)
logging.debug(str(num)+'calculated EMF3')

EMF_zeddy2        = EMF_zeddy_tmean + EMF_zeddy_teddy
logging.debug(str(num)+'calculated EMF4')


EMF = {'EMF_zeddy': EMF_zeddy, 'EMF_zeddy2': EMF_zeddy2, \
       'EMF_zeddy_tmean': EMF_zeddy_tmean, 'EMF_zeddy_teddy': EMF_zeddy_teddy}


def div_EMF(X):
  cos_phi      = np.deg2rad(np.cos(lat))[None,None,:]
  div =  (1/(a*cos_phi**2)) * (np.gradient (X*cos_phi**2, lat, axis=-1 ))
  return div

div_EMF_zeddy        = div_EMF( zonal_eddy_flux(a=ucomp, b=vcomp))
logging.debug(str(num)+'calculated EMF1')

div_EMF_zeddy_tmean  = div_EMF( zonal_eddy_time_mean_flux(a=ucomp, b=vcomp))
logging.debug(str(num)+'calculated EMF2')

div_EMF_zeddy_teddy  = div_EMF( zonal_eddy_time_eddy_flux(a=ucomp, b=vcomp))
logging.debug(str(num)+'calculated EMF3')

div_EMF_zeddy2        = div_EMF( EMF_zeddy_tmean + EMF_zeddy_teddy)
logging.debug(str(num)+'calculated EMF4')


divEMF = {'div_EMF_zeddy': div_EMF_zeddy, 'div_EMF_zeddy2': div_EMF_zeddy2, \
          'div_EMF_zeddy_tmean': div_EMF_zeddy_tmean, 'div_EMF_zeddy_teddy': div_EMF_zeddy_teddy}



u_w_zeddy_tmean =  zonal_eddy_time_mean_flux( a=ucomp, b=omega )
logging.debug(str(num)+'calculated uw1')

u_w_zeddy_teddy =  zonal_eddy_time_eddy_flux( a=ucomp, b=omega )
logging.debug(str(num)+'calculated uw2')

u_w_zeddy       =  zonal_eddy_flux( a=ucomp, b=omega )
logging.debug(str(num)+'calculated uw3')

u_w_zeddy2      =  u_w_zeddy_tmean + u_w_zeddy_teddy
logging.debug(str(num)+'calculated uw4')


u_w = {'u_w_zeddy_tmean': u_w_zeddy_tmean, 'u_w_zeddy_teddy': u_w_zeddy_teddy, \
       'u_w_zeddy': u_w_zeddy, 'u_w_zeddy2': u_w_zeddy2}

## ugrad_dp(a=ucomp, ax=1)

du_w_zeddy_tmean_dp =  ugrad_dp(a=u_w_zeddy_tmean, ax=1)
logging.debug(str(num)+'calculated duw1_dp')

du_w_zeddy_teddy_dp =  ugrad_dp(a=u_w_zeddy_teddy, ax=1)
logging.debug(str(num)+'calculated duw2_dp')

du_w_zeddy_dp       =  ugrad_dp(a=u_w_zeddy, ax=1)
logging.debug(str(num)+'calculated duw3_dp')

du_w_zeddy2_dp      =  ugrad_dp(a=u_w_zeddy2, ax=1)
logging.debug(str(num)+'calculated duw4_dp')


du_w_dp = {'du_w_zeddy_tmean_dp': du_w_zeddy_tmean_dp, 'du_w_zeddy_teddy_dp': du_w_zeddy_teddy_dp, \
           'du_w_zeddy_dp': du_w_zeddy_dp, 'du_w_zeddy2_dp': du_w_zeddy2_dp}

    
momentum_terms = {'fv_zmean_tmean':fv_zmean_tmean,'du_dt':du_dt, 'w_dudp':w_dudp, 'zeta_v':zeta_v, \
                  'EMF':EMF, 'divEMF':divEMF, 'u_w':u_w, 'du_w_dp':du_w_dp }


save(destination+"momentum_terms.hkl", momentum_terms)
save(destination+"coord.hkl", coord)

logging.debug(str(num)+"loaded momentum terms")
   
end = ti.time()

logging.debug(str(num)+"-Finished- time taken = "+str( (end-start))+" seconds" )  


