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

logging.basicConfig( filename = log_directory+'reload_momentum'+exp_dir+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
         
logging.debug("** INTERPOLATION ** ")
    
            
####################
#### soomthening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth            



destination = quarter+'/post_process_data/isca_repeat/avged_over'+str(no_of_days)+'days/'+exp_dir+num+'/'
make_sure_path_exists(destination);


v_MM_MMC            = []
v_TE_MMC            = []
v_SE_MMC            = []
v_net_MMC           = []

v_MM1_MMC            = []
v_TE1_MMC            = []
v_SE1_MMC            = []
v_net1_MMC           = []

v_MM_MMC_v            = []
v_TE_MMC_v            = []
v_SE_MMC_v            = []
v_net_MMC_v           = []

v_MM1_MMC_v            = []
v_TE1_MMC_v            = []
v_SE1_MMC_v            = []
v_net1_MMC_v           = []

v_MM_del_m_mmc       = []
v_SE_del_m_mmc       = []
v_TE_del_m_mmc       = []
v_net_del_m_mmc      = []

v_MM1_del_m_mmc       = []
v_SE1_del_m_mmc       = []
v_TE1_del_m_mmc       = []
v_net1_del_m_mmc      = []

v_MM_del_v_mmc       = []
v_SE_del_v_mmc       = []
v_TE_del_v_mmc       = []
v_net_del_v_mmc      = []

v_MM1_del_v_mmc       = []
v_SE1_del_v_mmc       = []
v_TE1_del_v_mmc       = []
v_net1_del_v_mmc      = []

v_MMC               = []
vel                 = []

v_tot = []
v_TE  = []
v_SE  = []
v_MM  = []
v_net = []
v_TE1  = []
v_SE1  = []
v_MM1  = []
v_net1 = []
zeta   = []

destination = quarter+'/post_process_data/isca_repeat/avged_over'+str(no_of_days)+'days/'+exp_dir+'/'
make_sure_path_exists(destination);


for i in range(0,int(num)):
  source=quarter+'/post_process_data/isca_repeat/avged_over'+str(no_of_days)+'days/'+exp_dir+str(i)+'/'
  
  v_MMC.append(load(source+"v_MMC.hkl"))
  vel.append(load(source+"vs.hkl"))
   
  v_MM_MMC.append(v_MMC[i]['v_MM_MMC']['v_MM_flux'])
  v_TE_MMC.append(v_MMC[i]['v_TE_MMC']['v_TE_flux'])
  v_SE_MMC.append(v_MMC[i]['v_SE_MMC']['v_SE_flux'])
  v_net_MMC.append(v_MMC[i]['v_net_MMC']['v_net_flux'])

  v_MM1_MMC.append(v_MMC[i]['v_MM1_MMC']['v_MM1_flux'])
  v_TE1_MMC.append(v_MMC[i]['v_TE1_MMC']['v_TE1_flux'])
  v_SE1_MMC.append(v_MMC[i]['v_SE1_MMC']['v_SE1_flux'])
  v_net1_MMC.append(v_MMC[i]['v_net1_MMC']['v_net1_flux'])

  v_MM_MMC_v.append(v_MMC[i]['v_MM_MMC']['v_MM_flux_vert'])
  v_TE_MMC_v.append(v_MMC[i]['v_TE_MMC']['v_TE_flux_vert'])
  v_SE_MMC_v.append(v_MMC[i]['v_SE_MMC']['v_SE_flux_vert'])
  v_net_MMC_v.append(v_MMC[i]['v_net_MMC']['v_net_flux_vert'])

  v_MM1_MMC_v.append(v_MMC[i]['v_MM1_MMC']['v_MM1_flux_vert'])
  v_TE1_MMC_v.append(v_MMC[i]['v_TE1_MMC']['v_TE1_flux_vert'])
  v_SE1_MMC_v.append(v_MMC[i]['v_SE1_MMC']['v_SE1_flux_vert'])
  v_net1_MMC_v.append(v_MMC[i]['v_net1_MMC']['v_net1_flux_vert'])
   
  v_MM_del_m_mmc.append(v_MMC[i]['v_MM_MMC']['v_MM_del_m_mmc'])
  v_TE_del_m_mmc.append(v_MMC[i]['v_TE_MMC']['v_TE_del_m_mmc'])
  v_SE_del_m_mmc.append(v_MMC[i]['v_SE_MMC']['v_SE_del_m_mmc'])
  v_net_del_m_mmc.append(v_MMC[i]['v_net_MMC']['v_net_del_m_mmc'])

  v_MM1_del_m_mmc.append(v_MMC[i]['v_MM1_MMC']['v_MM1_del_m_mmc'])
  v_TE1_del_m_mmc.append(v_MMC[i]['v_TE1_MMC']['v_TE1_del_m_mmc'])
  v_SE1_del_m_mmc.append(v_MMC[i]['v_SE1_MMC']['v_SE1_del_m_mmc'])
  v_net1_del_m_mmc.append(v_MMC[i]['v_net1_MMC']['v_net1_del_m_mmc'])

  v_MM_del_v_mmc.append(v_MMC[i]['v_MM_MMC']['v_MM_del_v_mmc'])
  v_TE_del_v_mmc.append(v_MMC[i]['v_TE_MMC']['v_TE_del_v_mmc'])
  v_SE_del_v_mmc.append(v_MMC[i]['v_SE_MMC']['v_SE_del_v_mmc'])
  v_net_del_v_mmc.append(v_MMC[i]['v_net_MMC']['v_net_del_v_mmc'])

  v_MM1_del_v_mmc.append(v_MMC[i]['v_MM1_MMC']['v_MM1_del_v_mmc'])
  v_TE1_del_v_mmc.append(v_MMC[i]['v_TE1_MMC']['v_TE1_del_v_mmc'])
  v_SE1_del_v_mmc.append(v_MMC[i]['v_SE1_MMC']['v_SE1_del_v_mmc'])
  v_net1_del_v_mmc.append(v_MMC[i]['v_net1_MMC']['v_net1_del_v_mmc'])


  v_tot.append(vel[i]['v_tot'])
  v_TE.append(vel[i]['v_TE'])
  v_SE.append(vel[i]['v_SE'])
  v_MM.append(vel[i]['v_MM'])
  v_net.append(vel[i]['v_net'])
  v_TE1.append(vel[i]['v_TE1'])
  v_SE1.append(vel[i]['v_SE1'])  
  v_MM1.append(vel[i]['v_MM1'])  
  v_net1.append(vel[i]['v_net1'])
  zeta.append(vel[i]['zeta'])
  
  logging.debug("----"+str(i)+"----")
 
coord=(load(source+"coord_dic.hkl"))
lat      =coord['lat']
lon      =coord['lon']
no_of_plevels=coord["no_of_plevels"]


def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS           = np.shape(v_tot)[1]
print MONTHS

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

def zon_int_vert(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None,None]))*a
    return y/10**15

import scipy.integrate as integrate
              

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


v_tot              = interpolate_vertical(M(v_tot));     
v_TE               = interpolate_vertical(M(v_TE));
v_SE               = interpolate_vertical(M(v_SE));
v_MM               = interpolate_vertical(M(v_MM));
v_net              = interpolate_vertical(M(v_net));
v_TE1               = interpolate_vertical(M(v_TE1));
v_SE1               = interpolate_vertical(M(v_SE1));
v_MM1               = interpolate_vertical(M(v_MM1));
v_net1              = interpolate_vertical(M(v_net1));
zeta                = interpolate_vertical(M(zeta));

v_MM_MMC_v        = interpolate_vertical(M(v_MM_MMC_v)); 
v_SE_MMC_v        = interpolate_vertical(M(v_SE_MMC_v));
v_TE_MMC_v        = interpolate_vertical(M(v_TE_MMC_v));       
v_net_MMC_v       = interpolate_vertical(M(v_net_MMC_v));

v_MM1_MMC_v        = interpolate_vertical(M(v_MM1_MMC_v));         
v_SE1_MMC_v        = interpolate_vertical(M(v_SE1_MMC_v));
v_TE1_MMC_v        = interpolate_vertical(M(v_TE1_MMC_v));   
v_net1_MMC_v       = interpolate_vertical(M(v_net1_MMC_v));

v_MM_MMC        = interpolate(M(v_MM_MMC));
v_SE_MMC        = interpolate(M(v_SE_MMC));
v_TE_MMC        = interpolate(M(v_TE_MMC));
v_net_MMC       = interpolate(M(v_net_MMC));

v_MM1_MMC        = interpolate(M(v_MM1_MMC));
v_SE1_MMC        = interpolate(M(v_SE1_MMC));
v_TE1_MMC        = interpolate(M(v_TE1_MMC));
v_net1_MMC       = interpolate(M(v_net1_MMC));

v_MM_del_m_mmc        = interpolate(M(v_MM_del_m_mmc));
v_SE_del_m_mmc        = interpolate(M(v_SE_del_m_mmc));
v_TE_del_m_mmc        = interpolate(M(v_TE_del_m_mmc));
v_net_del_m_mmc        = interpolate(M(v_net_del_m_mmc));

v_MM1_del_m_mmc        = interpolate(M(v_MM1_del_m_mmc));
v_SE1_del_m_mmc        = interpolate(M(v_SE1_del_m_mmc));
v_TE1_del_m_mmc        = interpolate(M(v_TE1_del_m_mmc));
v_net1_del_m_mmc        = interpolate(M(v_net1_del_m_mmc));

v_MM_del_v_mmc        = interpolate(M(v_MM_del_v_mmc));
v_SE_del_v_mmc        = interpolate(M(v_SE_del_v_mmc));
v_TE_del_v_mmc        = interpolate(M(v_TE_del_v_mmc));
v_net_del_v_mmc        = interpolate(M(v_net_del_v_mmc));

v_MM1_del_v_mmc        = interpolate(M(v_MM1_del_v_mmc));
v_SE1_del_v_mmc        = interpolate(M(v_SE1_del_v_mmc));
v_TE1_del_v_mmc        = interpolate(M(v_TE1_del_v_mmc));
v_net1_del_v_mmc        = interpolate(M(v_net1_del_v_mmc));


logging.debug("Calculated interpolation")

##################################################
#########      Save in dictionaries      ##########
###################################################

v_decomp   = {'v_tot':v_tot,'v_TE':v_TE, 'v_SE':v_SE, 'v_MM':v_MM, 'v_net':v_net,\
              'v_TE1':v_TE, 'v_SE1':v_SE1, 'v_MM1':v_MM1, 'v_net1':v_net1, 'zeta':zeta}


v_MMC      = {'v_MM_MMC':v_MM_MMC, 'v_SE_MMC': v_SE_MMC, 'v_TE_MMC': v_TE_MMC, \
              'v_MM1_MMC':v_MM1_MMC, 'v_SE1_MMC': v_SE1_MMC, 'v_TE1_MMC': v_TE1_MMC,\
              'v_MM_MMC_v':v_MM_MMC_v, 'v_SE_MMC_v': v_SE_MMC_v, 'v_TE_MMC_v': v_TE_MMC_v, \
              'v_MM1_MMC_v':v_MM1_MMC_v, 'v_SE1_MMC_v': v_SE1_MMC_v, 'v_TE1_MMC_v': v_TE1_MMC_v,\
              'v_net_MMC':v_net_MMC, 'v_net1_MMC':v_net1_MMC,\
              'v_net_MMC_v':v_net_MMC_v, 'v_net1_MMC_v':v_net1_MMC_v}


thermo_dynamic = {'v_MM_del_m_mmc':v_MM_del_m_mmc, 'v_MM_del_v_mmc':v_MM_del_v_mmc,\
                  'v_TE_del_m_mmc':v_TE_del_m_mmc, 'v_TE_del_v_mmc':v_TE_del_v_mmc,\
                  'v_SE_del_m_mmc':v_SE_del_m_mmc, 'v_SE_del_v_mmc':v_SE_del_v_mmc,\
                  'v_net_del_m_mmc':v_net_del_m_mmc, 'v_net_del_v_mmc':v_net_del_v_mmc,\
                  'v_MM1_del_m_mmc':v_MM1_del_m_mmc, 'v_MM1_del_v_mmc':v_MM1_del_v_mmc,\
                  'v_TE1_del_m_mmc':v_TE1_del_m_mmc, 'v_TE1_del_v_mmc':v_TE1_del_v_mmc,\
                  'v_SE1_del_m_mmc':v_SE1_del_m_mmc, 'v_SE1_del_v_mmc':v_SE1_del_v_mmc,\
		  'v_net1_del_m_mmc':v_net1_del_m_mmc, 'v_net1_del_v_mmc':v_net1_del_v_mmc}


save(destination+"v_decomp.hkl", v_decomp)
save(destination+"coord.hkl", coord)
save(destination+"v_MMC.hkl", v_MMC)

logging.debug("Saved averaged momentum terms")

end = ti.time()
logging.debug("Looks great !! Time taken --> "+str(end-start))


