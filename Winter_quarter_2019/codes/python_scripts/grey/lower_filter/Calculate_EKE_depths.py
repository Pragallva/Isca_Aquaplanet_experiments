import numpy as np
import numpy as np
from scipy.signal import freqz
from scipy.signal import butter, lfilter


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

filename=str(dirc[1])
num=str(dirc[2])

log_directory='/project2/tas1/pragallva/Winter_quarter_2019/codes/shell_script/greyT85/log/'+'EKE'+filename
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


#dirc=sys.argv

source='/project2/tas1/pragallva/Winter_quarter_2019/exp_data/grey/'+filename+'/'
one_year=source+filename+num+'.nc'
destination='/project2/tas1/pragallva/Winter_quarter_2019/post_process_data/grey/EKE/'+filename+num+'/'
make_sure_path_exists(destination);

v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('...........Calculate EKE.............'+str(num))

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

time_shift=int(-1440*(1-0.75))
def shift(x):
    return np.roll(x, time_shift, axis=0)


sigma_full=v_var['pfull'][::-1]
sigma_half=v_var['phalf'][::-1]
sigma_half[-1]=0.0001
p_sfc=shift(v_var['ps'][:]) ## Pa
lat=v_var['lat'][:]
lon=v_var['lon'][:]

time= nc.num2date(times, units=v_var['time'].units, calendar= v_var['time'].calendar )

logging.debug('imported coordinates')


no_of_years = len(time)/(12*30*4)
no_of_months =12
days = 30
hours = 4

pres_full=shift(v_var['pres_full'][:,::-1,:,:])
pres_half=shift(v_var['pres_half'][:,::-1,:,:])

logging.debug('imported pressure')


v_comp    =shift(v_var['vcomp'][:,::-1,:,:])
logging.debug('saved v_comp')
u_comp    =shift(v_var['ucomp'][:,::-1,:,:])
logging.debug('saved u_comp')
Z         =shift(v_var['height'][:,::-1,:,:])
logging.debug('saved Z')


#### Filtering signal ##
logging.debug('Filtering signal now')

#HOURS=60.0*60.0
#DAYS = 24 * HOURS
    # Sample rate and desired cutoff frequencies (in Hz).
#fs = 4.0/(24.0)
#lowcut =  4.0/(10*24.0)
#highcut = 4.0/(2*24.0) 


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


def FILTER(data):
    lowcut  = 4.0/(6.5*24) ## 1/days
    highcut = 4.0/(2.0*24) ## 1/days
    order=3
    b, a = butter(order, [lowcut,highcut], btype='band')
    filtered_signal = lfilter(b, a, data)
    return filtered_signal

u_prime=np.copy(u_comp)
v_prime=np.copy(v_comp)
Z_prime=np.copy(Z)

for p in range(len(sigma_full)):
    for la in range(len(lat)):
        for lo in range(len(lon)):
            u_prime[:,p,la,lo]= FILTER(u_comp[:,p,la,lo])
            v_prime[:,p,la,lo]= FILTER(v_comp[:,p,la,lo])
            Z_prime[:,p,la,lo]= FILTER(Z[:,p,la,lo]-Z[:,p,la,lo].mean())
logging.debug('Filtering u and v at plevel '+str(sigma_full[p]/1000))


u_sq = u_prime**2
v_sq = v_prime**2
Z_sq = Z_prime**2
EKE  = (u_sq+v_sq)/2.0

EMF = u_prime*v_prime

logging.debug('calculated EKE')


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


vert_EKE  = vertical_integral(EKE )
vert_EMF  =  vertical_integral(EMF)
vert_Zsq  = vertical_integral(Z_sq)


def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(sigma_full),len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)

EKE_dic      ={'v_sq':R(v_sq),'u_sq':R(u_sq),'Z_sq':R(Z_sq),'vert_Zsq':vert_Zsq,'EKE':R(EKE), 'vert_EKE':vert_EKE, 'EMF': R(EMF), 'vert_EMF':vert_EMF, 'lat':lat, 'lon':lon, 'sigma_full': sigma_full}

logging.debug('saved EKE data')
save(destination+"EKE.hkl",EKE_dic)
logging.debug('EKE calculated !! Great !!')