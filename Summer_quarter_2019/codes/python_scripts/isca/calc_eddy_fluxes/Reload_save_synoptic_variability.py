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
FREQ       = dirc[3]
days       = int(dirc[4])
no_of_days = days

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

logging.basicConfig( filename = log_directory+'synoptic_interp'+exp_dir+FREQ+str(days)+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
         
logging.debug("** INTERPOLATION ** ")
    
            
####################
#### soomthening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth            


EKE                 = []
vert_EKE            = []

EMF                 = []
vert_EMF            = []

EKE_flux            = []
sensible_flux       = []
latent_flux         = []
potential_flux      = []

vert_EKE_flux       = []
vert_sensible_flux  = []
vert_latent_flux    = []
vert_potential_flux = []

momentum_fluxes     = []

destination = quarter+'/post_process_data/isca_repeat/avged_over'+str(no_of_days)+'days/'+exp_dir+'/'
make_sure_path_exists(destination);

for i in range(0,int(num)):
  source = quarter+'/post_process_data/isca_repeat/avged_over'+str(no_of_days)+'days/'+exp_dir+str(i)+'/'
  
  momentum_fluxes.append(load(source+"eddies_"+FREQ+'_freq_'+'avg'+str(days)+'.hkl'))

  EKE.append(momentum_fluxes[i]['EKE'])
  vert_EKE.append(momentum_fluxes[i]['vert_EKE'])
 
  EMF.append(momentum_fluxes[i]['EMF'])
  vert_EMF.append(momentum_fluxes[i]['vert_EMF'])
 
  EKE_flux.append(momentum_fluxes[i]['EKE_flux'])
  sensible_flux.append(momentum_fluxes[i]['sensible_flux'])
  latent_flux.append(momentum_fluxes[i]['latent_flux'])
  potential_flux.append(momentum_fluxes[i]['potential_flux'])

  vert_sensible_flux.append(momentum_fluxes[i]['vert_sensible_flux'])
  vert_latent_flux.append(momentum_fluxes[i]['vert_latent_flux'])
  vert_potential_flux.append(momentum_fluxes[i]['vert_potential_flux'])
  vert_EKE_flux.append(momentum_fluxes[i]['vert_EKE_flux'])

  logging.debug("----"+str(i)+"----")
 

coord        = momentum_fluxes[0]
lat          = coord['lat']
lon          = coord['lon']
pres         = coord['sigma_full']

def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

LAT=len(lat)
a=6371.0e3

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS           = np.shape(EKE_flux)[1]

if int(MONTHS)==360/days :
   logging.debug("Month calculation looks good")

print MONTHS

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

def zon_int_vert(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None,None]))*a
    return y/10**15

import scipy.integrate as integrate
              
no_of_plevels = len(pres)       

def interpolate_vertical(X):
    interp=np.zeros((LATN,no_of_plevels,MONTHS))
    for m in range(MONTHS):
        for p in range(no_of_plevels):
            interpolation_function = interp1d(lat, X[m,p,:],kind='linear')
            interp[:,p,m]=interpolation_function(latn)
    return interp
  
def interpolate(X):
     interp=np.zeros((LATN,MONTHS))
     for m in range(MONTHS):
             interpolation_function = interp1d(lat, X[m,:],kind='linear')
             interp[:,m]=interpolation_function(latn)
     return interp


       
EKE                 = interpolate_vertical(M(EKE))
vert_EKE            = interpolate(M(vert_EKE))

EMF                 = interpolate_vertical(M(EMF))
vert_EMF            = interpolate(M(vert_EMF))

EKE_flux            = interpolate_vertical(M(EKE_flux))
sensible_flux       = interpolate_vertical(M(sensible_flux))
latent_flux         = interpolate_vertical(M(latent_flux))
potential_flux      = interpolate_vertical(M(potential_flux))

vert_sensible_flux  = interpolate(M(vert_sensible_flux))
vert_latent_flux    = interpolate(M(vert_latent_flux))
vert_potential_flux = interpolate(M(vert_potential_flux))
vert_EKE_flux       = interpolate(M(vert_EKE_flux))

logging.debug("Calculated interpolation")

##################################################
#########      Save in dictionaries      ##########
###################################################

def R(x):
   return x


eddies_dic       = {'EKE':R(EKE), 'vert_EKE':vert_EKE, 'EMF': R(EMF), 'vert_EMF':vert_EMF,'EKE_flux':R(EKE_flux), \
                   'sensible_flux':R(sensible_flux),          'latent_flux':R(latent_flux),        'potential_flux': R(potential_flux), \
                   'vert_sensible_flux':vert_sensible_flux,         'vert_latent_flux':vert_latent_flux, 'vert_potential_flux':vert_potential_flux, \
                   'vert_EKE_flux':R(vert_EKE_flux),'lat':latn, 'lon':lon, 'sigma_full': pres, 'frequency': FREQ, 'avg_days': days}

save(destination+"eddies_"+FREQ+'_freq_'+'avg'+str(days)+'.hkl', eddies_dic)

logging.debug("Saved averaged synoptic eddies")

end = ti.time()
logging.debug("Looks great !! Time taken --> "+str(end-start))


