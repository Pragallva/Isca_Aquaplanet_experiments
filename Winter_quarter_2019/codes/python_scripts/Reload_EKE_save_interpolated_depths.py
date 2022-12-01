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

exp=str(dirc[1])
name=str(dirc[2])
num=str(dirc[3]) ## represents an annual year of data

log_directory='/project2/tas1/pragallva/Winter_quarter_2019/codes/shell_script/log/'+name
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

# represents an annual year of data
#source='/project2/tas1/pragallvaring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'

destination='/project2/tas1/pragallva/Winter_quarter_2019/post_process_data/'+name+'/'
make_sure_path_exists(destination);

EKE=[]
raw_data=[]
EKE_vert=[]
EMF=[]
EMF_vert=[]
Z_sq=[]
Z_sq_vert=[]


#EKE_dic      ={'EKE':R(EKE), 'vert_EKE':vert_EKE, 'lat':lat, 'lon':lon, 'sigma_full': sigma_full}


for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Winter_quarter_2019/post_process_data/'+name+str(i)+'/'
 raw_data.append(load(source+"EKE.hkl"))    
 EKE.append(raw_data[i]['EKE'])
 EKE_vert.append(raw_data[i]['vert_EKE'])
 EMF.append(raw_data[i]['EMF'])
 EMF_vert.append(raw_data[i]['vert_EMF'])
 Z_sq.append(raw_data[i]['Z_sq'])
 Z_sq_vert.append(raw_data[i]['vert_Zsq'])
    
 logging.debug("----"+str(i)+"----")

coord_source = '/project2/tas1/pragallva/Fall_quarter_2019/post_process_data/HC10_la50m_oc5m0/' 
coord=(load(coord_source+"coord_dic.hkl"))
lat      =coord['lat']
lon      =coord['lon']
no_of_plevels=coord["no_of_plevels"]

def M(x):
    return np.mean(x,axis=0)

EKE=M(EKE) ; EKE_vert=M(EKE_vert)
EMF=M(EMF) ; EMF_vert=M(EMF_vert)
Z_sq=M(Z_sq); Z_sq_vert=M(Z_sq_vert)

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
EMF_vert_interp =interpolate(EMF_vert)
EMF_interp      =interpolate_plevels(EMF)
Zsq_interp      =interpolate_plevels(Z_sq)
Zsq_vert_interp =interpolate(Z_sq_vert)

logging.debug("Calculated interpolation")

raw_data_dict={'Z_sq':Zsq_interp,'Zsq_vert':Zsq_vert_interp,'EKE':EKE_interp,'EKE_vert':EKE_vert_interp, 'EMF':EMF_interp,'EMF_vert':EMF_vert_interp, 'latn':latn,'sigma_full':raw_data[0]['sigma_full'],'latn':latn}
save(destination+'EKE_interp.hkl', raw_data_dict)
logging.debug("Saved raw data interp dictionary")

# import shutil
# for i in range(0,int(num)):
#  source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
#  shutil.rmtree(source)

# logging.debug("Deleted individual 30 years !!")

logging.debug("Looks great !!")


