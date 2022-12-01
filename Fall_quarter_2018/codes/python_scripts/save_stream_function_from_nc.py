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

logging.debug('imported nc file')

lat=v_var['lat'][:]
lon=v_var['lon'][:]

logging.debug('imported coordinates')

# In[10]:

v_comp    =v_var['vcomp'][:,::-1,:,:]
logging.debug('opened v_comp again')
u_comp    =v_var['ucomp'][:,::-1,:,:]
logging.debug('opened u_comp again')

from windspharm.standard import VectorWind
from windspharm.tools    import prep_data, recover_data, order_latdim
from windspharm.examples import example_data_path

# Read zonal and meridional wind components from file using the netCDF4
# module. The components are defined on pressure levels and are in separate
# files.

uwnd = u_comp
vwnd = v_comp

lons = lon
lats = lat

# The standard interface requires that latitude and longitude be the leading
# dimensions of the input wind components, and that wind components must be
# either 2D or 3D arrays. The data read in is 3D and has latitude and
# longitude as the last dimensions. The bundled tools can make the process of
# re-shaping the data a lot easier to manage.
uwnd, uwnd_info = prep_data(uwnd, 'tyx')
vwnd, vwnd_info = prep_data(vwnd, 'tyx')

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
sf, vp = w.sfvp()

sf = recover_data(sf, uwnd_info)
vp = recover_data(vp, uwnd_info)
logging.debug('Create a recover data')

stream = {"stream_func":sf,"velocity_pot":vp}

save(destination+"stream_dic.hkl", stream)
logging.debug("saved stream function dictionary")





