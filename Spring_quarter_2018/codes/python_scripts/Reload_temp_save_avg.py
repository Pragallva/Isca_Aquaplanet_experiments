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
log_directory='/project2/tas1/pragallva/Spring_quarter_2018/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
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
#source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'

destination='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'
make_sure_path_exists(destination);

TEMP=[]
temp=[]

for i in range(1,int(num)+1):
 source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
 TEMP.append(load(source+"temp_dic.hkl"))
 temp.append(TEMP[i-1]['temp'])
 logging.debug("----"+str(i)+"----")

def M(x):
    return np.mean(x,axis=0)

temp=M(temp) 
temp_dic={'temp':temp.mean(axis=2).mean(axis=2).mean(axis=0)}

save(destination+'temp_dic.hkl', temp_dic)

logging.debug("Looks great !!")


