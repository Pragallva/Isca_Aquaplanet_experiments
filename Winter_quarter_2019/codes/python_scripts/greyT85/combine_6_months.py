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
import glob
import logging

filename  = str(dirc[1])
num       = str(dirc[2])
pres_surf = int(dirc[3])

log_directory='/project2/tas1/pragallva/Winter_quarter_2019/codes/shell_script/greyT85/log/'+filename+'_combine'
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

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
            
####################
#### soomthening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth    


durations=['first6months/*hkl','last6months/*hkl']

if pres_surf==1 :
    source='/project2/tas1/pragallva/Winter_quarter_2019/post_process_data/greyT85/data_in_pres_coord/'+filename+str(num)+'/'
else:
    source='/project2/tas1/pragallva/Winter_quarter_2019/post_process_data/greyT85/'+filename+str(num)+'/'

dict1 = (glob.glob(source+durations[0]))
dict2 = (glob.glob(source+durations[1]))   

def combine(d1,d2,fieldname):
      final={}
      for key in d1.keys():
         if isinstance(d1[key], (list, tuple, np.ndarray)) :
            temporary={key: np.concatenate( (d1[key], d2[key]), axis=0) }
            final.update(temporary)
            logging.debug("updated  "+str(key))
         else :
            temporary2={key: (d1[key])}  
            final.update(temporary2)
            logging.debug("updated  "+str(key))
                  
      save(source+str(fieldname), final)
      logging.debug("..... combined and saved "+str(fieldname)+' ......')

for i in range(len(dict1)):
    fieldname = os.path.split(dict1[i])[1]
    D1 = dict1[i]; D2= dict2[i]
    d1=load(D1); d2=load(D2)
    if fieldname=='coord_dic.hkl':    
        save(source+str(fieldname), d1)
    else: 
        combine(d1,d2,fieldname)
    
       
logging.debug("Looks great !!")


