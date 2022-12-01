from netCDF4 import Dataset, num2date
import metpy.calc as mcalc
from metpy.cbook import get_test_data
from metpy.plots import add_metpy_logo
from metpy.units import units
import numpy as np
import netCDF4 as nc
import sys
import os.path
import pylab as py
import matplotlib.cm as cm
import numpy.ma as ma

import sys
import os
import errno
dirc=sys.argv

import logging


import time as ti
start = ti.time()

exp_dir  = str(dirc[1]) ## "HC0_la50m_oc50m"
num      = str(dirc[2]) ## represents an annual year of data
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

logging.basicConfig( filename = log_directory+exp_dir+'stress.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

source   = quarter+'/exp_data/isca_repeat/'+exp_dir+'/'
one_year = source+exp_dir+num+'.nc'


destination = quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+num+'/'
make_sure_path_exists(destination);

data = nc.Dataset(one_year,'r')

logging.debug(str(num)+'...........extract_surface_stress.py.............')

logging.debug(str(num)+'imported nc file')

def get_data(variable='transient', f = one_year):
    ncfile = f
    v_var  = nc.Dataset(ncfile,'r')
    data   = v_var.variables
    each_model_name = data[variable][:]
    v_var.close()
    return each_model_name

tau_u       = get_data('tau_u')
logging.debug(str(num)+'imported tau u')

pres_sur    = get_data('ps')
logging.debug(str(num)+'imported ps')

pres_full = get_data('pres_full') 
logging.debug(str(num)+'imported pressure full')

pres_diff1 = pres_sur - pres_full[:,5,:,:]
pres_diff2 = pres_sur - pres_full[:,1,:,:]

g = 10.0
dtau_by_pres_diff1 = -g*tau_u/(pres_diff1)
dtau_by_pres_diff2 = -g*tau_u/(pres_diff2)
dtau_by_surf_pres  = -g*tau_u/(pres_sur)


lat      = get_data('lat')
lon      = get_data('lon')
time     = get_data('time')


if no_of_days == None:
   no_of_days = 30
else:
   days = int(no_of_days)

no_of_months = 360/days
hours        = 4
no_of_years  = len(time)/(no_of_months*days*hours)


def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)

def reshape(y):
    y1 = y.reshape((no_of_years,no_of_months,days,hours,len(lat),len(lon)))
    return y1


def zonal_mean(X):
    return np.nanmean(X, axis=-1)[..., None]

def time_mean(X):
    print X.shape
    daily_mean   = np.nanmean(X, axis=3)[:,:,:,None,:,:]          ## Averaged across all hours
    monthly_mean = np.nanmean(daily_mean, axis=2)[:,:,None,:,:,:] ## Averaged across all months
    return monthly_mean

def time_eddy(X):
    return X-time_mean(X)

def zonal_eddy(X):
    return X-time_mean(X)


def combo(X):
    return zonal_mean(time_mean(reshape(X)))

logging.debug(str(num)+'imported coordinates')

### Create a stress file ###

stress={'tau_u': combo(tau_u),'ps': combo(pres_sur),'dtau_by_pres_diff1':combo(dtau_by_pres_diff1),'dtau_by_pres_diff2': combo(dtau_by_pres_diff2),\
        'dtau_by_surf_pres':combo(dtau_by_surf_pres),'lat': lat,'lon':lon }

save(destination+'surf_stress.hkl',stress)

logging.debug(str(num)+'Saved the stress file')   

end = ti.time()

logging.debug(str(num)+'Program ends')
logging.debug(str(num)+'Total time taken - '+str(end-start))
logging.debug(str(num)+'------------------------')
