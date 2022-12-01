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
log_directory='/project2/tas1/pragallva/Summer_quarter_2018/codes/shell_script/log/'+'SFCtemp'+dirc[1]+'_'+dirc[2]
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


#dirc=sys.argv
num=str(dirc[3]) ## represents an annual year of data
source='/project2/tas1/pragallva/Summer_quarter_2018/exp_data/'+dirc[1]+'_'+dirc[2]+'/'
one_year=source+dirc[1]+'_'+dirc[2]+num+'.nc'

destination='/project2/tas1/pragallva/Summer_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+num+'/'
make_sure_path_exists(destination);
v_variables = nc.Dataset(one_year,'r')
v_var=v_variables.variables

logging.debug('...........Extracting surface temperature.............'+str(num))

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
sigma_full=v_var['pfull'][::-1]
sigma_half=v_var['phalf'][::-1]
sigma_half[-1]=0.0001
p_sfc=v_var['ps'][:]  ## Pa
lat=v_var['lat'][:]
lon=v_var['lon'][:]

time= nc.num2date(times, units=v_var['time'].units, calendar= v_var['time'].calendar )

logging.debug('imported coordinates')


no_of_years = len(time)/(12*30*4)
no_of_months =12
days = 30
hours = 4

t_surf = v_var['t_surf'][:,:,:]

logging.debug('imported t_surf')

def R(y):
    y1=y.reshape((no_of_years,no_of_months,days,hours,len(lat),len(lon)))
    return y1.mean(axis=0).mean(axis=1).mean(axis=1).mean(axis=-1)
                  

tsurf_dic      ={'tsurf':R(t_surf), 'lat':lat, 'lon':lon}


logging.debug('saved tsurf data')

save(destination+"tsurf.hkl",tsurf_dic)

logging.debug('t_surf saved !! Great !!')
