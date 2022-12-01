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

def read_small_field(field_name):
    field = data.variables[field_name][:]
    return field

def read_big_field(field_name):
    UNITS = units(data.variables[field_name].units)
    if field_name=='temp':
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
     
    num=str(0) ## represents an annual year of data
    dirc=np.array([0,'aqua','isca5m'])
    source='/project2/tas1/pragallva/Spring_quarter_2018/exp_data/'+dirc[1]+'_'+dirc[2]+'/'
    one_year=source+dirc[1]+'_'+dirc[2]+num+'.nc'
    data = nc.Dataset(one_year,'r')
    
    ########## Read data sets ############
    field_names_without_unit=['lat','lon', 'pfull', 'time','ucomp','vcomp', 'height', 'temp', 'sphum', 'pres_full']
    #field_names_with_unit=['ucomp','vcomp', 'height', 'temp', 'sphum', 'pres_full']
    
    for name in field_names_without_unit:
        y=read_small_field(name)
        print str(name), y.shape

    #############################################





