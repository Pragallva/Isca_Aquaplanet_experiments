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
from multiprocessing import Pool, current_process, cpu_count, Process, Array, Queue

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

def read_small_field(field_name,field):    
    for index, fname in enumerate(field_name):
        field.put(data.variables[fname][:])
        print fname
    
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
    
    field_names=['lat','lon','ucomp']
    #, 'pfull', 'time','ucomp','vcomp', 'height', 'temp', 'sphum', 'pres_full']
    
    num=str(0) ## represents an annual year of data
    dirc=np.array([0,'aqua','isca5m'])
    source='/project2/tas1/pragallva/Spring_quarter_2018/exp_data/'+dirc[1]+'_'+dirc[2]+'/'
    one_year=source+dirc[1]+'_'+dirc[2]+num+'.nc'
    data = nc.Dataset(one_year,'r')
    
    field=Queue()
    p=Process(target=read_small_field,args=(field_names,field))
    
    p.start()
    p.join()
    
    while field.empty() is False:
        print(field.get().shape)
    
    
    
    ########## Read data sets ############
    field_names_without_unit=['lat','lon', 'pfull', 'time','ucomp','vcomp', 'height', 'temp', 'sphum', 'pres_full']
    
    
    
"""
    pool = Pool()
    lat  = pool.map(read_small_field,('lat'))
    lon  = pool.map(read_small_field,('lon'))
    u    = pool.map(read_small_field,('ucomp'))
    v    = pool.map(read_small_field,('vcomp'))
    ht   = pool.map(read_small_field,('height'))

#    u= u.get(timeout=10)
#    v= v.get(timeout=10)
#    ht= ht.get(timeout=10)
    print u.shape, v.shape, ht.shape
    
    pool.close()
    pool.join()
    
"""
        
"""    
    pool_size = len(field_names_without_unit)#+len(field_names_with_unit) # log this value
    pool_read = Pool(processes = pool_size, initializer=start_process)
    lat = pool_read.map( read_small_field, field_names_without_unit )
#    u, v, ht, temp, sphum, pres = pool_read.map( read_big_field, field_names_with_unit )
  # , lon, sigma, times, u, v, ht, temp, sphum, pres
    pool_read.close()
    pool_read.join()
    
    print np.shape. v.shape. temp.shape, lat.shape, pool_size
    #############################################
"""




