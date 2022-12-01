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
#log_directory='/project2/tas1/pragallva/Summer_quarter_2018/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
#logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

import hickle as hkl

from multiprocessing import Pool, current_process, cpu_count

def start_process():
    print("start:",current_process().name)

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

def logging_fun():
    return

def read_small_field(field_name):
    field = data.variables[field_name][:]
    return field

def read_big_field(field_name):
    UNITS = units(data.variables[field_name].units)
    if field_name='temp':
        UNITS=units.K
    field = data.variables[field_name][:]  *  UNITS
    return field
    
def isobar_levels(field,arg):
    field = ma.masked_invalid(isobaric_levels[arg].magnitude)
    return field

def reshape(x): # average over years
    x=x.reshape(no_of_years,no_of_months,days,hours,len(plevs),len(lat),len(lon))
    return x.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)


if __name__=='__main__':
     
    num=str(dirc[3]) ## represents an annual year of data
    source='/project2/tas1/pragallva/Spring_quarter_2018/exp_data/'+dirc[1]+'_'+dirc[2]+'/'
    one_year=source+dirc[1]+'_'+dirc[2]+num+'.nc'

    destination='/project2/tas1/pragallva/Summer_quarter_2018/post_process_data/data_in_pres_coord/'+dirc[1]+'_'+dirc[2]+num+'/'
    make_sure_path_exists(destination);
    data = nc.Dataset(one_year,'r')
    
    ########## Read data sets ############
    field_names_without_unit=['lat','lon', 'pfull', 'time']
    field_names_with_unit=['ucomp','vcomp', 'height', 'temp', 'sphum', 'pres_full']
    
    pool_size=len(field_names_without_unit)+len(field_names_with_unit) # log this value
    pool_read = Pool(processes = pool_size, initializer=start_process)
    lat, lon, sigma, times = pool_read.map( read_small_field, field_names_without_unit )
    u, v, ht, temp, sphum, pres = pool_read.map( read_big_field, field_names_with_unit )
    pool_read.close()
    pool_read.join()
    
    print u.shape. v.shape. temp.shape, lat.shape
    #############################################
    

#logging.debug(str(num)+'...........interpolation_from_sigma_to_pres_coord.py.............')

#logging.debug(str(num)+'imported nc file')

lat   = data.variables['lat'][:]
lon   = data.variables['lon'][:]
sigma  = data.variables['pfull'][:] # sigma_coord
height= data.variables['height'][:]
times =data.variables['time'][:]

logging.debug(str(num)+'imported coordinates, height')

u = data.variables['ucomp'][:]   * units(data.variables['ucomp'].units)
logging.debug(str(num)+'imported u')
v = data.variables['vcomp'][:]   * units(data.variables['vcomp'].units)
logging.debug(str(num)+'imported v')
ht = data.variables['height'][:] * units(data.variables['height'].units)
logging.debug(str(num)+'imported ht')
temp = data.variables['temp'][:] * units.K
logging.debug(str(num)+'imported temp')
sphum = data.variables['sphum'][:] * units(data.variables['sphum'].units)
logging.debug(str(num)+'imported sphum')
pres = data.variables['pres_full'][:] * units(data.variables['pres_full'].units)
logging.debug(str(num)+'imported pressure')

time= nc.num2date(times, units=data.variables['time'].units, calendar= data.variables['time'].calendar )
# goes one time to main function

no_of_years = len(time)/(12*30*4)
no_of_months =12
days = 30
hours = 4

# Array of desired pressure levels
plevs = [0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0] * units.hPa

isobaric_levels = mcalc.log_interp(plevs, pres, ht, temp, u , v , sphum, axis=1)
## goes to main function

logging.debug(str(num)+'interpolation done')


def aoy(x): # average over years
    #x=x.reshape(1,x.shape[0]/1,x.shape[1],x.shape[2],x.shape[3])
    return x#.mean(axis=0)

def isobar_levels(field,arg):
    field = ma.masked_invalid(isobaric_levels[arg].magnitude)
    return field

htp    = aoy( ma.masked_invalid(isobaric_levels[0].magnitude) )
temp   = aoy( ma.masked_invalid(isobaric_levels[1].magnitude) )
up     = aoy( ma.masked_invalid(isobaric_levels[2].magnitude) )
vp     = aoy( ma.masked_invalid(isobaric_levels[3].magnitude) )
sphump = aoy( ma.masked_invalid(isobaric_levels[4].magnitude) )

# Save zonally averaged #

no_of_years = len(time)/(12*30*4)
no_of_months =12
days = 30
hours = 4
def reshape(x): # average over years
    x=x.reshape(no_of_years,no_of_months,days,hours,len(plevs),len(lat),len(lon))
    return x.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)

vel    = {"u" :reshape(up), "v":reshape(vp),      "plevs":plevs.magnitude}
others = {"ht":reshape(htp),"temp":reshape(temp), "sphum": reshape(sphump)}

logging.debug(str(num)+'formed dictionaries')

save(destination+"u_v.hkl", vel)
save(destination+"ht_temp_sphum.hkl", others)

logging.debug(str(num)+"saved velocity and others in pressure coordinates")

"""
from windspharm.standard import VectorWind
from windspharm.tools    import prep_data, recover_data, order_latdim
from windspharm.examples import example_data_path

# Read zonal and meridional wind components from file using the netCDF4
# module. The components are defined on pressure levels and are in separate
# files.
import numpy.ma as ma

uwnd = up.filled(fill_value=0)
vwnd = vp.filled(fill_value=0)

lons = lon
lats = lat

# The standard interface requires that latitude and longitude be the leading
# dimensions of the input wind components, and that wind components must be
# either 2D or 3D arrays. The data read in is 3D and has latitude and
# longitude as the last dimensions. The bundled tools can make the process of
# re-shaping the data a lot easier to manage.
uwnd, uwnd_info = prep_data(uwnd, 'tzyx') ## Note that tzyx indicates time, level, lat, lon
vwnd, vwnd_info = prep_data(vwnd, 'tzyx')

logging.debug('reordering uwind and vwind')

# It is also required that the latitude dimension is north-to-south. Again the
# bundled tools make this easy.
lats, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)

# Create a VectorWind instance to handle the computation of streamfunction and
# velocity potential.
w = VectorWind(uwnd, vwnd)
logging.debug('Create a VectorWind instance')

# Compute the streamfunction and velocity potential. Also use the bundled
# tools to re-shape the outputs to the 4D shape of the wind components as they
# were read off files.
sf, velp = w.sfvp()

sf   = recover_data(sf, uwnd_info)
velp = recover_data(velp, uwnd_info)

logging.debug('VectorWind recovered')

sf=sf.squeeze()[:,:,::-1,:]
### This windspharm library needs data from North to south. But my data is from south to north, hence I am reversing here.
sf= ( ma.masked_array(sf, up.mask) )
velp=velp.squeeze()[:,:,::-1,:]
velp= ( ma.masked_array(velp, up.mask) )

logging.debug('Masked the data')

stream = {"stream_func":sf,"velocity_pot":velp}

save(destination+"stream_dic.hkl", stream)
"""




