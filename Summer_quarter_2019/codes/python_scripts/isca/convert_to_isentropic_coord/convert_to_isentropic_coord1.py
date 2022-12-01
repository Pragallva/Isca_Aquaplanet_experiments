from netCDF4 import Dataset, num2date
import metpy.calc as mcalc
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

def as_str(x):
   for k, v in list(locals().iteritems()):
       if v is x:
             x_as_str = k
   return x_as_str

log_directory = quarter+'/codes/shell_scripts/log/'
make_sure_path_exists( path   = log_directory )

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig( filename = log_directory+exp_dir+'_isen1.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

source   = quarter+'/exp_data/isca_repeat/data_in_pres_coord/'+exp_dir+'/'
one_year = source+exp_dir+num+'.nc'

destination = quarter+'/exp_data/isca_repeat/data_in_isentropic_coord/'+exp_dir+'/'
make_sure_path_exists(destination);

data = nc.Dataset(one_year,'r')

logging.debug(str(num)+'...........interpolation_from_pres_to_isentrope_coord.py.............')

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

plevs    = get_data('level')*100 * units.pascals

lat      = get_data('lat')
lon      = get_data('lon')
time     = get_data('time')

logging.debug(str(num)+'imported coordinates')

# Array of desired pressure levels
isen_levels = np.arange(270,420,5) * units.K

def interp_to_isentrope(field1=height, field2=temp): 
        t=0
	field1_isen = np.zeros((len(time),len(isen_levels),field1.shape[2],field1.shape[3]))
        field2_isen = np.copy(field1_isen)

        isen1_levels1 = mcalc.isentropic_interpolation(isen_levels, plevs, temp[...], \
                                                  field1[...], field2[...], axis=1)
	logging.debug(str(num)+'interpolation done')

	field1_isen[...]     =  np.squeeze(isen1_levels1[0].magnitude)
	field2_isen[...]     =  np.squeeze(isen1_levels1[1].magnitude)
   
        return np.squeeze(field1_isen), np.squeeze(field2_isen) 


logging.debug('Calculating for height and temperature isentrope')
hp, tp = interp_to_isentrope( height, temp)
logging.debug('Calculation successful: shape of hp and tp '+str(np.shape(hp)) )

#qp, up = interp_to_isentrope( sphum, ucomp)

#vp, wp = interp_to_isentrope( vcomp, omega)

### Create a netcdf file ###

## create a dimension
dataset   = nc.Dataset(destination+exp_dir+num+'test2.nc','w',format='NETCDF4_CLASSIC')
levels    = dataset.createDimension('level', len(isen_levels))
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
#ucomps  = dataset.createVariable('ucomp',  np.float32, ('time','level','lat','lon'))
#vcomps  = dataset.createVariable('vcomp',  np.float32, ('time','level','lat','lon'))
#omegas  = dataset.createVariable('omega',  np.float32, ('time','level','lat','lon'))
#sphums  = dataset.createVariable('sphum',  np.float32, ('time','level','lat','lon'))
heights = dataset.createVariable('height', np.float32, ('time','level','lat','lon'))
logging.debug(str(num)+'Create variables')

lats[:]    = lat
lons[:]    = lon
times[:]   = time
levels[:]  = isen_levels.magnitude

temps[:,:,:,:]   = tp
#ucomps[:,:,:,:]  = up
#vcomps[:,:,:,:]  = vp
#omegas[:,:,:,:]  = wp
#sphums[:,:,:,:]  = qp
heights[:,:,:,:] = hp

logging.debug(str(num)+'assigned dataset')

dataset.close()

end = ti.time()

logging.debug(str(num)+'Program ends')
logging.debug(str(num)+'Total time taken - '+str(end-start))
logging.debug(str(num)+'------------------------')
