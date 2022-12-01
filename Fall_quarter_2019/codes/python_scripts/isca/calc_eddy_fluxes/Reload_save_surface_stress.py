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
import time as ti

start = ti.time()

exp_dir    = str(dirc[1]) ## "HC0_la50m_oc50m"
num        = str(dirc[2]) ## represents an annual year of data
no_of_days = dirc[3]


current_dir = os.getcwd()
parts       = current_dir.split(os.sep)
quarter     = "/"+os.path.join(*parts[:5])


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


log_directory = quarter+'/codes/shell_scripts/log/'
make_sure_path_exists( path   = log_directory )

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig( filename = log_directory+'eddy_interp'+exp_dir+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
         
logging.debug("** INTERPOLATION ** ")
    
            
####################
#### soomthening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth            

tau_u               = []
ps                  = []
dtau_by_pres_diff1  = []
dtau_by_pres_diff2  = []
dtau_by_surf_pres   = []


momentum_fluxes     = []


destination = quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+'/'
make_sure_path_exists(destination);


for i in range(0,int(num)):
  source=quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+str(i)+'/'
  
  momentum_fluxes.append(load(source+"surf_stress.hkl"))

  tau_u.append(np.squeeze(momentum_fluxes[i]['tau_u']))
  ps.append(np.squeeze(momentum_fluxes[i]['ps']))
  dtau_by_pres_diff1.append(np.squeeze(momentum_fluxes[i]['dtau_by_pres_diff1']))
  dtau_by_pres_diff2.append(np.squeeze(momentum_fluxes[i]['dtau_by_pres_diff2']))
  dtau_by_surf_pres.append(np.squeeze(momentum_fluxes[i]['dtau_by_surf_pres']))

  logging.debug("----"+str(i)+"----")
 

coord        =(load(source+"coord.hkl"))
lat          = coord['lat']
lon          = coord['lon']

def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS           = np.shape(dtau_by_pres_diff1)[1]
print MONTHS

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

def zon_int_vert(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None,None]))*a
    return y/10**15

import scipy.integrate as integrate
              

def interpolate_vertical(X):
    interp=np.zeros((LATN,MONTHS))
    for m in range(MONTHS):
            interpolation_function = interp1d(lat, X[m,:],kind='linear')
            interp[:,m]=interpolation_function(latn)
    return interp
         

tau_u                  = interpolate_vertical(M(tau_u));     
ps                     = interpolate_vertical(M(ps));
dtau_by_pres_diff1     = interpolate_vertical(M(dtau_by_pres_diff1));         
dtau_by_pres_diff2     = interpolate_vertical(M(dtau_by_pres_diff2));
dtau_by_surf_pres      = interpolate_vertical(M(dtau_by_surf_pres));

logging.debug("Calculated interpolation")

##################################################
#########      Save in dictionaries      ##########
###################################################


surf_stress = {'tau_u': tau_u, 'ps': ps, \
               'dtau_by_pres_diff1': dtau_by_pres_diff1, \
               'dtau_by_pres_diff2': dtau_by_pres_diff2, \
               'dtau_by_surf_pres': dtau_by_surf_pres}


save(destination+"surf_stress.hkl", surf_stress)

logging.debug("Saved averaged surface stress terms")

end = ti.time()
logging.debug("Looks great !! Time taken --> "+str(end-start))


