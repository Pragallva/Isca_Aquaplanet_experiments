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

edge=int(dirc[1])
depth=int(dirc[2])
num=int(dirc[3]) ## represents an annual year of data

log_directory='/project2/tas1/pragallva/Winter_quarter_2019/codes/shell_script/log/'+'HC'+str(edge)+'_la'+str(depth)+'m_oc'+str(depth)+'m'
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

destination='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/data_in_pres_coord/'+'HC%d_la%dm_oc%dm/'%(edge,depth,depth)
make_sure_path_exists(destination);

raw_data = []

psi_prime_weighted      = []
psi_prime_weighted_v = []

psi_prime_unweighted      = [] 
psi_prime_unweighted_v = []

pot_prime_weighted      = []
pot_prime_weighted_v = []

pot_prime_unweighted      = []
pot_prime_unweighted_v    = []


for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/data_in_pres_coord/'+'HC%d_la%dm_oc%dm%d/'%(edge,depth,depth,i)
 raw_data.append(load(source+"stream_rms.hkl"))    
 psi_prime_weighted.append(raw_data[i]['psi_prime_weighted'])
 psi_prime_weighted_v.append(raw_data[i]['psi_prime_weighted_v'])
 psi_prime_unweighted.append(raw_data[i]['psi_prime_unweighted'])
 psi_prime_unweighted_v.append(raw_data[i]['psi_prime_unweighted_v'])  
 pot_prime_weighted.append(raw_data[i]['pot_prime_weighted'])
 pot_prime_weighted_v.append(raw_data[i]['pot_prime_weighted_v'])
 pot_prime_unweighted.append(raw_data[i]['pot_prime_unweighted'])
 pot_prime_unweighted_v.append(raw_data[i]['pot_prime_unweighted_v'])

 logging.debug("----"+str(i)+"----")

coord_source = '/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/HC10_la50m_oc5m0/' 
coord=(load(coord_source+"coord_dic.hkl"))
lat      =coord['lat']
lon      =coord['lon']

def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

psi_prime_weighted=M(psi_prime_weighted) ; psi_prime_weighted_v=M(psi_prime_weighted_v)
psi_prime_unweighted=M(psi_prime_unweighted) ; psi_prime_unweighted_v=M(psi_prime_unweighted_v)

pot_prime_weighted=M(pot_prime_weighted) ; pot_prime_weighted_v=M(pot_prime_weighted_v)
pot_prime_unweighted=M(pot_prime_unweighted) ; pot_prime_unweighted_v=M(pot_prime_unweighted_v)

no_of_plevels = psi_prime_weighted_v.shape[1]

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS=12
        
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



psi_prime_weighted =interpolate(psi_prime_weighted)
psi_prime_unweighted = interpolate(psi_prime_unweighted)
pot_prime_weighted =interpolate(pot_prime_weighted)
pot_prime_unweighted = interpolate(pot_prime_unweighted)

psi_prime_weighted_v =interpolate_plevels(psi_prime_weighted_v)
psi_prime_unweighted_v = interpolate_plevels(psi_prime_unweighted_v)
pot_prime_weighted_v =interpolate_plevels(pot_prime_weighted_v)
pot_prime_unweighted_v = interpolate_plevels(pot_prime_unweighted_v)

logging.debug("Calculated interpolation")

stream_rms = {"psi_prime_weighted":psi_prime_weighted,"psi_prime_weighted_v":psi_prime_weighted_v, "psi_prime_unweighted":psi_prime_unweighted, "psi_prime_unweighted_v":psi_prime_unweighted_v,"pot_prime_weighted":pot_prime_weighted, "pot_prime_weighted_v":pot_prime_weighted_v, "pot_prime_unweighted": pot_prime_unweighted, "pot_prime_unweighted_v":pot_prime_unweighted_v, 'latn':latn }

save(destination+'stream_rms.hkl', stream_rms)
logging.debug("Saved interp stream_rms dictionary")

# import shutil
# for i in range(0,int(num)):
#  source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
#  shutil.rmtree(source)

# logging.debug("Deleted individual 30 years !!")

logging.debug("Looks great !!")


