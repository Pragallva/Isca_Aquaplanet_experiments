import numpy as np
import numpy as np
from scipy.signal import freqz
from scipy.signal import butter, lfilter, filtfilt
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
FREQ       = str(dirc[3])
days       = int(dirc[4])
ORDER      = 3

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

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


log_directory = quarter+'/codes/shell_scripts/log/'
make_sure_path_exists( path   = log_directory )

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig( filename = log_directory+'synoptic'+exp_dir+FREQ+str(days)+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

source   = quarter+'/exp_data/am2_repeat/'+exp_dir+'/'
one_year = source+exp_dir+num+'.nc'
pk_bk_source = quarter+'/exp_data/am2_repeat/pk_bk_source.nc'

destination = quarter+'/post_process_data/am2_repeat/'+exp_dir+num+'/'
make_sure_path_exists(destination);


v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('...........calculate synoptic variability.............'+str(num))

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


def get_data(variable='ucomp', f = one_year, T=1):
    ncfile = f
    v_var  = nc.Dataset(ncfile,'r')
    data   = v_var.variables
    each_model_name = data[variable][:]
    v_var.close()
    if T==1:
        return (each_model_name.transpose(2,1,0,3))[:1440,::-1,:,:]
    else:
        return each_model_name


ucomp    =get_data('ucomp') #* units.meter / (units.second)
logging.debug(str(num)+'imported u')

vcomp    = get_data('vcomp') # * units.meter / (units.second)
logging.debug(str(num)+'imported v')

Z       = get_data('zfull')  #* units.meter
logging.debug(str(num)+'imported Z')

q      = get_data('sphum') #* units.dimensionless
logging.debug(str(num)+'imported q')

temp     = get_data('temp')# * units.K
logging.debug(str(num)+'imported temp')

p_sfc  = get_data(variable='ps',T=0).transpose(1,0,2)[:1440,:,:]
logging.debug(str(num)+'imported surface pressure')



lat      = get_data(variable='lat',T=0)
lon      = get_data(variable='lon',T=0)
time     = get_data(variable='time',T=0)[:1440]

sigma_full = get_data(variable='pfull',T=0)[::-1]/1000.0
sigma_half = get_data(variable='phalf',T=0)[::-1]/1000.0

pk       = get_data('pk',pk_bk_source,0)[::-1]
bk       = get_data('bk',pk_bk_source,0)[::-1]

logging.debug(str(num)+'imported coordinates')


def convert_sigma_to_pressure_coordinates():
    ## stack surface pressure to all 26 pressure levels
    ps=np.stack([p_sfc]*len(bk), axis=1)
    phalf=(pk[None,:,None,None])+(bk[None,:,None,None])*ps ## Because sigma came with 1000 already multiplied to it
    pfull=(phalf[:,:-1,:,:]-phalf[:,1:,:,:])/(np.log(phalf[:,:-1,:,:])-np.log(phalf[:,1:,:,:]))
    return phalf,pfull

pres_half=convert_sigma_to_pressure_coordinates()[0] #* units.pascals
logging.debug(str(num)+'phalf calculated')
pres_full=convert_sigma_to_pressure_coordinates()[1] #* units.pascals
logging.debug(str(num)+'pfull calculated')


logging.debug('imported coordinates')

if days == None :
   days = 30
no_of_years = len(time)/(12*30*4)
no_of_months = 360/days
hours = 4

logging.debug('imported pressure')


m = Cp*temp + g*Z + L*q

#### Filtering signal ##
logging.debug('Filtering signal now')


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y




def FILTER(data,low_day=2.0, high_day=10.0, band_type='high_pass'):
    lowcut  = 4.0/(high_day*24) ## 1/days
    highcut = 4.0/(low_day*24) ## 1/days
    order   = ORDER
    if ((band_type == 'high_pass') or (band_type == 'monthly_pass')):
        b, a = butter(order, lowcut, btype='high')
    else :
        b, a = butter(order, [lowcut, highcut], btype='band')
    filtered_signal = filtfilt(b, a, data)
    return filtered_signal


def filterall(Y=ucomp, freq='high'):
    if (freq == 'high_pass'):
        low_day = 2.0;
        high_day = 10.0;
    if (freq == 'monthly_pass'):
        low_day = 2.0;
        high_day = 30.0;
    if (freq == 'high'):
        low_day = 2.0;
        high_day = 10.0
    if (freq == 'inter1'):
        low_day  = 10.0;
        high_day = 30.0;
    if (freq == 'inter2'):
        low_day  = 30.0;
        high_day = 45.0;
    if (freq == 'low'):
        low_day  = 45.0;
        high_day = 90.0;
    prime      = np.copy(Y)
    for p in range(len(sigma_full)):
      for la in range(len(lat)):
        for lo in range(len(lon)):
            prime[:,p,la,lo]= FILTER(Y[:,p,la,lo], low_day, high_day, band_type=freq)
      logging.debug('Filtering at'+str(sigma_full[p])+' sigma level for '+freq+' frequency')
    return prime




logging.debug('Filtering U')
u_prime   = filterall(ucomp, FREQ)

logging.debug('Filtering V')
v_prime   = filterall(vcomp, FREQ)

logging.debug('Filtering CpT')
CpT_prime = filterall(Cp*temp, FREQ)

logging.debug('Filtering gZ')
gZ_prime  = filterall(g*Z, FREQ)

logging.debug('Filtering Lq')
Lq_prime  = filterall(L*q, FREQ)

u_sq      = u_prime**2
v_sq      = v_prime**2

EKE  = (u_sq+v_sq)/2.0
EMF  = u_prime*v_prime

sensible_flux_prime  = CpT_prime*v_prime
latent_flux_prime    = Lq_prime*v_prime
potential_flux_prime = gZ_prime*v_prime 
EKE_flux_prime       = EKE*v_prime
logging.debug('calculated all transient fluxes')
                  
                  
def reshape_pres(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_half),len(lat),len(lon)))
    return y1
                  
def reshape(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1

ps=reshape_pres(pres_half)
weights = (ps[:,:,:,:,:-1,:,:]-ps[:,:,:,:,1:,:,:])/g

def weighted(arg):
    w = arg*weights
    return w

def vertical_integral(m):
    M         = weighted(reshape(m))                 # (year, month,days, hour, plev, lat, lon)        
    monthly_M = M.mean(axis=2).mean(axis=2)          # (year, month, plev, lat, lon)
    zonal_M   = monthly_M.mean(axis=-1)              # (year, month, plev, lat)
    vert_M    = zonal_M.sum(axis=2)                  # (year, month, lat)
    return vert_M.mean(axis=0)                       # (month, lat)--> yearly average


vert_EKE            = vertical_integral(EKE)
vert_EMF            = vertical_integral(EMF)
vert_latent_flux    = vertical_integral(latent_flux_prime)
vert_sensible_flux  = vertical_integral(sensible_flux_prime)
vert_potential_flux = vertical_integral(potential_flux_prime)
vert_EKE_flux       = vertical_integral(EKE_flux_prime)

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return weighted(y1).mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)

eddies_dic      ={'EKE':R(EKE), 'vert_EKE':vert_EKE, 'EMF': R(EMF), 'vert_EMF':vert_EMF,'EKE_flux':R(EKE_flux_prime),\
                  'sensible_flux':R(sensible_flux_prime), 'latent_flux':R(latent_flux_prime), 'potential_flux': R(potential_flux_prime), \
                  'vert_sensible_flux':vert_sensible_flux, 'vert_latent_flux':vert_latent_flux, 'vert_potential_flux':vert_potential_flux, \
                  'vert_EKE_flux':vert_EKE_flux,'lat':lat, 'lon':lon, 'sigma_full': sigma_full, 'frequency': FREQ, 'avg_days': days}

logging.debug('saved eddies data')
save(destination+"eddies_"+FREQ+'_freq_'+'avg'+str(days)+'.hkl',eddies_dic)
end = ti.time()
logging.debug('Eddies calculated !! Great !!')
logging.debug('Time taken = '+str(end-start))


