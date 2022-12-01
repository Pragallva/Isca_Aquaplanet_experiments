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
log_directory='/project2/tas1/pragallva/Fall_quarter_2017/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

###############################
# try:                        #
#    import cPickle as pickle #       Unfortunately pickle doesn't work
# except:                     #
#    import pickle            #
###############################

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
source='/project2/tas1/pragallva/Fall_quarter_2017/exp_data/'+dirc[1]+'/'
one_year=source+dirc[1]+'_'+dirc[2]+'.nc'

destination='/project2/tas1/pragallva/Fall_quarter_2017/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'
make_sure_path_exists(destination);
v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables


py.rc('text', usetex=True)
py.rc('font', family='serif', serif='Palatino',weight='bold')

filename= one_year #"/project2/tas1/pragallva/Fall_quarter_2017/exp_data/land/land_real.nc"
land_file="/project2/tas1/pragallva/Fall_quarter_2017/land_files/Ruth_full.nc"

data=nc.Dataset(filename,'r')
land_data=nc.Dataset(land_file,'r')

logging.debug('...........interpolation_from_sigma_to_pres_coord.py.............')

logging.debug('imported nc file')

lat   = data.variables['lat'][:]
lon   = data.variables['lon'][:]
sigma  = data.variables['pfull'][:] # sigma_coord
height= data.variables['height'][:]

logging.debug('imported coordinates, height')

u = data.variables['ucomp'][:]   * units(data.variables['ucomp'].units)
logging.debug('imported u')
v = data.variables['vcomp'][:]   * units(data.variables['vcomp'].units)
logging.debug('imported v')
ht = data.variables['height'][:] * units(data.variables['height'].units)
logging.debug('imported ht')
temp = data.variables['temp'][:] * units.K
logging.debug('imported temp')
sphum = data.variables['sphum'][:] * units(data.variables['sphum'].units)
logging.debug('imported sphum')
pres = data.variables['pres_full'][:] * units(data.variables['pres_full'].units)
logging.debug('imported pressure')

# Array of desired pressure levels
plevs = [0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0] * units.hPa

isobaric_levels = mcalc.log_interp(plevs, pres, ht, temp, u , v , sphum, axis=1)

logging.debug('interpolation done')

def aoy(x): # average over years
    x=x.reshape(5,x.shape[0]/5,x.shape[1],x.shape[2],x.shape[3])
    return x.mean(axis=0)

htp    = aoy( ma.masked_invalid(isobaric_levels[0].magnitude) )
temp   = aoy( ma.masked_invalid(isobaric_levels[1].magnitude) )
up     = aoy( ma.masked_invalid(isobaric_levels[2].magnitude) )
vp     = aoy( ma.masked_invalid(isobaric_levels[3].magnitude) )
sphump = aoy( ma.masked_invalid(isobaric_levels[4].magnitude) )

logging.debug('internal extraction after interpolation')


vel    = {"u":up,"v":vp}
others = {"ht":htp,"temp":temp, "sphum": sphump}
save(destination+"u_v.hkl", vel)
save(destination+"ht_temp_sphum.hkl", others)

logging.debug("saved raw coordinates dictionary")

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

logging.debug("saved stream function dictionary")


