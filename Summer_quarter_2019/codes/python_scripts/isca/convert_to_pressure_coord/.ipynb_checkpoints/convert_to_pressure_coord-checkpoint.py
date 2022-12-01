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

logging.basicConfig( filename = log_directory+exp_dir+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

source   = quarter+'/exp_data/isca_repeat/'+exp_dir+'/'
one_year = source+exp_dir+num+'.nc'

destination = quarter+'/exp_data/isca_repeat/data_in_pres_coord/'+exp_dir+'/'
make_sure_path_exists(destination);

data = nc.Dataset(one_year,'r')

logging.debug(str(num)+'...........interpolation_from_sigma_to_pres_coord.py.............')

logging.debug(str(num)+'imported nc file')

def get_data(variable='transient', f = one_year):
    ncfile = f
    v_var  = nc.Dataset(ncfile,'r')
    data   = v_var.variables
    each_model_name = data[variable][:]
    v_var.close()
    return each_model_name

ucomp    = get_data('ucomp')* units.meter / (units.second) 
logging.debug(str(num)+'imported u')

vcomp    = get_data('vcomp')* units.meter / (units.second) 
logging.debug(str(num)+'imported v')

omega    = get_data('omega')* units.pascals / (units.second)
logging.debug(str(num)+'imported w')

height   = get_data('height')* units.meter
logging.debug(str(num)+'imported Z')

sphum    = get_data('sphum')* units.dimensionless
logging.debug(str(num)+'imported q')

temp     = get_data('temp')* units.K
logging.debug(str(num)+'imported temp')

pres_full = get_data('pres_full') * units.pascals
logging.debug(str(num)+'imported pressure full')

pres_half = get_data('pres_half') * units.pascals
logging.debug(str(num)+'imported pressure half')


lat      = get_data('lat')
lon      = get_data('lon')
sigma    = np.append(get_data('pfull'), 1000.0)
time     = get_data('time')

logging.debug(str(num)+'imported coordinates')

#time= nc.num2date(times, units=data.variables['time'].units, calendar= data.variables['time'].calendar )

# Array of desired pressure levels
plevs = plevs = [0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0,  700.0, 750.0, 800.0, 825.0, 850.0, 900.0, 925.0, 950.0, 975.0, 1000.0] * units.hPa

isobaric_levels1 = mcalc.log_interp(plevs, pres_full, height, temp, axis=1)
logging.debug(str(num)+'interpolation done ht, T')

isobaric_levels2 = mcalc.log_interp(plevs, pres_full, sphum, ucomp, axis=1)
logging.debug(str(num)+'interpolation done q, U')

isobaric_levels3 = mcalc.log_interp(plevs, pres_full, vcomp, omega, axis=1)
logging.debug(str(num)+'interpolation done v, w')

hp     =  (isobaric_levels1[0].magnitude) 
logging.debug(str(num)+'height isobar')
tp     =  (isobaric_levels1[1].magnitude)
logging.debug(str(num)+'temp  isobar') 
qp     =  (isobaric_levels2[0].magnitude) 
logging.debug(str(num)+'moist  isobar') 
up     =  (isobaric_levels2[1].magnitude)
logging.debug(str(num)+'u  isobar')
vp     =  (isobaric_levels3[0].magnitude)
logging.debug(str(num)+'v  isobar')
wp     =  (isobaric_levels3[1].magnitude)
logging.debug(str(num)+'w  isobar')
### Create a netcdf file ###

## create a dimension
dataset   = nc.Dataset(destination+exp_dir+num+'.nc','w',format='NETCDF4_CLASSIC')
levels    = dataset.createDimension('level', len(plevs))
lats      = dataset.createDimension('lat', len(lat))
lons      = dataset.createDimension('lon', len(lon))
times     = dataset.createDimension('time', None)
logging.debug(str(num)+'Create dimension')   

## Netcdf attribute
dataset.title = exp_dir+"-"+num+" year"
lats     = dataset.createVariable('lat', np.float32, ('lat',))
lons     = dataset.createVariable('lon', np.float32, ('lon',))
levels   = dataset.createVariable('level', np.float32, ('level',))
times    = dataset.createVariable('time', np.float64, ('time',)) 
logging.debug(str(num)+'Create coordinates')

# Create the actual 4-d variable
    
temps   = dataset.createVariable('temp',   np.float32, ('time','level','lat','lon'))
ucomps  = dataset.createVariable('ucomp',  np.float32, ('time','level','lat','lon'))
vcomps  = dataset.createVariable('vcomp',  np.float32, ('time','level','lat','lon'))
omegas  = dataset.createVariable('omega',  np.float32, ('time','level','lat','lon'))
sphums  = dataset.createVariable('sphum',  np.float32, ('time','level','lat','lon'))
heights = dataset.createVariable('height', np.float32, ('time','level','lat','lon'))
logging.debug(str(num)+'Create variables')

lats[:]    = lat
lons[:]    = lon
times[:]   = time
levels[:]  = plevs.magnitude

temps[:,:,:,:]   = tp
ucomps[:,:,:,:]  = up
vcomps[:,:,:,:]  = vp
omegas[:,:,:,:]  = wp
sphums[:,:,:,:]  = qp
heights[:,:,:,:] = hp

logging.debug(str(num)+'assigned dataset')

dataset.close()

end = ti.time()

logging.debug(str(num)+'Program ends')
logging.debug(str(num)+'Total time taken - '+str(end-start))
logging.debug(str(num)+'------------------------')
