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
log_directory='/project2/tas1/pragallva/Summer_quarter_2018/codes/shell_script/log/SFC'+dirc[1]+'_'+dirc[2]
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

# try:
#    import cPickle as pickle        Unfortunately pickle doesn't work
# except:
#    import pickle

logging.debug("** INTERPOLATION ** ")

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
            
####################
#### soomthening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth            


num=str(dirc[3]) ## represents an annual year of data
#source='/project2/tas1/pragallvaring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'

destination='/project2/tas1/pragallva/Summer_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'
make_sure_path_exists(destination);

tsurf_data=[]
tsurf=[]

for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Summer_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
 tsurf_data.append(load(source+"tsurf.hkl"))    
 tsurf.append(tsurf_data[i]['tsurf'])
 
 logging.debug("----"+str(i)+"----")
 

def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

lat=tsurf_data[0]['lat']
lon=tsurf_data[0]['lon']
tsurf=M(tsurf) ; 

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS=12
tsurf_interp=np.zeros((LATN,MONTHS))

def interpolate(X):
    interp=np.zeros((LATN,MONTHS))
    for m in range(12):
            interpolation_function = interp1d(lat, X[m,:],kind='linear')
            interp[:,m]=interpolation_function(latn)
    return interp

tsurf_interp =interpolate(tsurf)

logging.debug("Calculated interpolation")


tsurf_dict={'tsurf':tsurf_interp,'latn':latn}
save(destination+'tsurf_interp.hkl', tsurf_dict)
logging.debug("Saved tsurf interp dictionary")

# import shutil
# for i in range(0,int(num)):
#  source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
#  shutil.rmtree(source)

# logging.debug("Deleted individual 30 years !!")

logging.debug("Looks great !!")


