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
log_directory='/project2/tas1/pragallva/Summer_quarter_2018/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
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

EKE=[]
raw_data=[]
EKE_vert=[]

#EKE_dic      ={'EKE':R(EKE), 'vert_EKE':vert_EKE, 'lat':lat, 'lon':lon, 'sigma_full': sigma_full}


for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Summer_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
 raw_data.append(load(source+"EKE.hkl"))    
 EKE.append(raw_data[i]['EKE'])
 EKE_vert.append(raw_data[i]['vert_EKE'])
    
 logging.debug("----"+str(i)+"----")
 
coord=(load(source+"coord_dic.hkl"))
lat      =coord['lat']
lon      =coord['lon']
no_of_plevels=coord["no_of_plevels"]

def M(x):
    return np.mean(x,axis=0)

EKE=M(EKE) ; EKE_vert=M(EKE_vert)

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS=12
EKE_interp=np.zeros((LATN,MONTHS))
EKE_vert_interp=np.zeros((LATN,no_of_plevels,MONTHS))

        
def interpolate_plevels(X):
    interp=np.zeros((LATN,no_of_plevels,MONTHS))
    for m in range(12):
        for p in range(no_of_plevels):
            interpolation_function = interp1d(lat, X[m,p,:],kind='linear')
            interp[:,p,m]=interpolation_function(latn)
    return interp

def interpolate(X):
    interp=np.zeros((LATN,MONTHS))
    for m in range(12):
            interpolation_function = interp1d(lat, X[m,:],kind='linear')
            interp[:,m]=interpolation_function(latn)
    return interp

EKE_vert_interp =interpolate(EKE_vert)
EKE_interp      =interpolate_plevels(EKE)


logging.debug("Calculated interpolation")


raw_data_dict={'EKE':EKE_interp,'EKE_vert':EKE_vert_interp,'latn':latn,'sigma_full':raw_data[0]['sigma_full'],'latn':latn}
save(destination+'EKE_interp.hkl', raw_data_dict)
logging.debug("Saved raw data interp dictionary")

# import shutil
# for i in range(0,int(num)):
#  source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
#  shutil.rmtree(source)

# logging.debug("Deleted individual 30 years !!")

logging.debug("Looks great !!")


