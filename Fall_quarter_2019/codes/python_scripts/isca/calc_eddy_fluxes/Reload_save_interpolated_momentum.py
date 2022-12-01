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



destination = quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+num+'/'
make_sure_path_exists(destination);



fv_zmean_tmean      = []
du_dt               = []

w_dudp_zmean_tmean  = []
w_dudp_zmean_teddy  = []
w_dudp_zmean        = []
w_dudp_zmean2       = []

zeta_v_zmean_tmean  = []
zeta_v_zmean_teddy  = []
zeta_v_zmean        = []
zeta_v_zmean2       = []
zeta_mean           = []

EMF_zeddy           = []
EMF_zeddy2          = []
EMF_zeddy_tmean     = []
EMF_zeddy_teddy     = []
EMF_zmean_teddy     = []
EMF_teddy           = []

div_EMF_zeddy       = []
div_EMF_zeddy2      = []
div_EMF_zeddy_tmean = []
div_EMF_zeddy_teddy = []
div_EMF_zmean_teddy = []
div_EMF_teddy       = []

u_w_zeddy_tmean     = []
u_w_zeddy_teddy     = []
u_w_zeddy           = []
u_w_zeddy2          = []
u_w_teddy           = []
u_w_zmean_teddy     = []

du_w_zeddy_tmean_dp = []
du_w_zeddy_teddy_dp = []
du_w_zeddy_dp       = []
du_w_zeddy2_dp      = []
du_w_teddy_dp       = []
du_w_zmean_teddy_dp = []

momentum_fluxes     = []
raw_data            = []

MSE  = []
CpT  = []
Lq   = []
gZ   = []
mv   = []
v    = []

destination = quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+'/'
make_sure_path_exists(destination);


for i in range(0,int(num)):
  source=quarter+'/post_process_data/isca_repeat/data_in_pres_coord/avged_over'+str(no_of_days)+'days/'+exp_dir+str(i)+'/'
  
  momentum_fluxes.append(load(source+"momentum_terms.hkl"))
  raw_data.append(load(source+"raw_data.hkl"))

  MSE.append(raw_data[i]['MSE'])
  CpT.append(raw_data[i]['CpT'])
  Lq.append(raw_data[i]['Lq'])
  gZ.append(raw_data[i]['gZ'])
  mv.append(raw_data[i]['mv'])
  dp_by_g = raw_data[i]['dp_by_g']  
  v.append(raw_data[i]['v'])  

  fv_zmean_tmean.append(momentum_fluxes[i]['fv_zmean_tmean'])
  du_dt.append(momentum_fluxes[i]['du_dt'])
 
  w_dudp_zmean_tmean.append(momentum_fluxes[i]['w_dudp']['w_dudp_zmean_tmean'])
  w_dudp_zmean_teddy.append(momentum_fluxes[i]['w_dudp']['w_dudp_zmean_teddy'])
  w_dudp_zmean.append(momentum_fluxes[i]['w_dudp']['w_dudp_zmean'])
  w_dudp_zmean2.append(momentum_fluxes[i]['w_dudp']['w_dudp_zmean2'])
 
  zeta_v_zmean_tmean.append(momentum_fluxes[i]['zeta_v']['zeta_v_zmean_tmean'])
  zeta_v_zmean_teddy.append(momentum_fluxes[i]['zeta_v']['zeta_v_zmean_teddy'])
  zeta_v_zmean.append(momentum_fluxes[i]['zeta_v']['zeta_v_zmean'])
  zeta_v_zmean2.append(momentum_fluxes[i]['zeta_v']['zeta_v_zmean2'])
  zeta_mean.append(momentum_fluxes[i]['zeta_v']['zeta_mean'])

  EMF_zeddy.append(momentum_fluxes[i]['EMF']['EMF_zeddy'])
  EMF_zeddy2.append(momentum_fluxes[i]['EMF']['EMF_zeddy2'])
  EMF_zeddy_tmean.append(momentum_fluxes[i]['EMF']['EMF_zeddy_tmean'])
  EMF_zeddy_teddy.append(momentum_fluxes[i]['EMF']['EMF_zeddy_teddy'])
  EMF_zmean_teddy.append(momentum_fluxes[i]['EMF']['EMF_zmean_teddy'])
  EMF_teddy.append(momentum_fluxes[i]['EMF']['EMF_teddy'])
   

  div_EMF_zeddy.append(momentum_fluxes[i]['divEMF']['div_EMF_zeddy'])
  div_EMF_zeddy2.append(momentum_fluxes[i]['divEMF']['div_EMF_zeddy2'])
  div_EMF_zeddy_tmean.append(momentum_fluxes[i]['divEMF']['div_EMF_zeddy_tmean'])
  div_EMF_zeddy_teddy.append(momentum_fluxes[i]['divEMF']['div_EMF_zeddy_teddy'])
  div_EMF_zmean_teddy.append(momentum_fluxes[i]['divEMF']['div_EMF_zmean_teddy'])
  div_EMF_teddy.append(momentum_fluxes[i]['divEMF']['div_EMF_teddy'])
  

  u_w_zeddy_tmean.append(momentum_fluxes[i]['u_w']['u_w_zeddy_tmean'])
  u_w_zeddy_teddy.append(momentum_fluxes[i]['u_w']['u_w_zeddy_teddy'])
  u_w_zeddy.append(momentum_fluxes[i]['u_w']['u_w_zeddy'])
  u_w_zeddy2.append(momentum_fluxes[i]['u_w']['u_w_zeddy2'])
  u_w_zmean_teddy.append(momentum_fluxes[i]['u_w']['u_w_zmean_teddy'])
  u_w_teddy.append(momentum_fluxes[i]['u_w']['u_w_teddy'])

  du_w_zeddy_tmean_dp.append(momentum_fluxes[i]['du_w_dp']['du_w_zeddy_tmean_dp'])
  du_w_zeddy_teddy_dp.append(momentum_fluxes[i]['du_w_dp']['du_w_zeddy_teddy_dp'])
  du_w_zeddy_dp.append(momentum_fluxes[i]['du_w_dp']['du_w_zeddy_dp'])
  du_w_zeddy2_dp.append(momentum_fluxes[i]['du_w_dp']['du_w_zeddy2_dp'])
  du_w_zmean_teddy_dp.append(momentum_fluxes[i]['du_w_dp']['du_w_zmean_teddy_dp'])
  du_w_teddy_dp.append(momentum_fluxes[i]['du_w_dp']['du_w_teddy_dp'])

  logging.debug("----"+str(i)+"----")
 

coord        =(load(source+"coord.hkl"))
lat          = coord['lat']
lon          = coord['lon']
pres         = coord['pres']

def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS           = np.shape(du_w_zeddy_tmean_dp)[1]
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
         

fv_zmean_tmean     = interpolate_vertical(M(fv_zmean_tmean));     
du_dt              = interpolate_vertical(M(du_dt));

w_dudp_zmean_tmean = interpolate_vertical(M(w_dudp_zmean_tmean)); 
w_dudp_zmean_teddy = interpolate_vertical(M(w_dudp_zmean_teddy));
w_dudp_zmean       = interpolate_vertical(M(w_dudp_zmean));       
w_dudp_zmean2      = interpolate_vertical(M(w_dudp_zmean2));

zeta_v_zmean_tmean = interpolate_vertical(M(zeta_v_zmean_tmean)); 
zeta_v_zmean_teddy = interpolate_vertical(M(zeta_v_zmean_teddy));
zeta_v_zmean       = interpolate_vertical(M(zeta_v_zmean));
zeta_v_zmean2      = interpolate_vertical(M(zeta_v_zmean2));
zeta_mean          = interpolate_vertical(M(zeta_mean));

EMF_zeddy          = interpolate_vertical(M(EMF_zeddy));          
EMF_zeddy2         = interpolate_vertical(M(EMF_zeddy2));
EMF_zeddy_teddy    = interpolate_vertical(M(EMF_zeddy_teddy));    
EMF_zeddy_tmean    = interpolate_vertical(M(EMF_zeddy_tmean));
EMF_zmean_teddy    = interpolate_vertical(M(EMF_zmean_teddy));
EMF_teddy          = interpolate_vertical(M(EMF_teddy));


div_EMF_zeddy       = interpolate_vertical(M(div_EMF_zeddy));       
div_EMF_zeddy2      = interpolate_vertical(M(div_EMF_zeddy2));
div_EMF_zeddy_tmean = interpolate_vertical(M(div_EMF_zeddy_tmean));  
div_EMF_zeddy_teddy = interpolate_vertical(M(div_EMF_zeddy_teddy));
div_EMF_zmean_teddy = interpolate_vertical(M(div_EMF_zmean_teddy));
div_EMF_teddy       = interpolate_vertical(M(div_EMF_teddy));

u_w_zeddy_tmean     = interpolate_vertical(M(u_w_zeddy_tmean));      
u_w_zeddy_teddy     = interpolate_vertical(M(u_w_zeddy_teddy));
u_w_zeddy           = interpolate_vertical(M(u_w_zeddy));           
u_w_zeddy2          = interpolate_vertical(M(u_w_zeddy2));
u_w_zmean_teddy     = interpolate_vertical(M(u_w_zmean_teddy));
u_w_teddy           = interpolate_vertical(M(u_w_teddy));

du_w_zeddy_tmean_dp = interpolate_vertical(M(du_w_zeddy_tmean_dp));  
du_w_zeddy_teddy_dp = interpolate_vertical(M(du_w_zeddy_teddy_dp));
du_w_zeddy_dp       = interpolate_vertical(M(du_w_zeddy_dp));        
du_w_zeddy2_dp      = interpolate_vertical(M(du_w_zeddy2_dp));
du_w_zmean_teddy_dp = interpolate_vertical(M(du_w_zmean_teddy_dp));
du_w_teddy_dp       = interpolate_vertical(M(du_w_teddy_dp));

MSE                 = interpolate_vertical(M(MSE));
CpT                 = interpolate_vertical(M(CpT));
Lq                  = interpolate_vertical(M(Lq));
gZ                  = interpolate_vertical(M(gZ));
mv                  = interpolate_vertical(M(mv));
v                   = interpolate_vertical(M(v));

logging.debug("Calculated interpolation")

##################################################
#########      Save in dictionaries      ##########
###################################################


u_w = {'u_w_zeddy_tmean': u_w_zeddy_tmean, 'u_w_zeddy_teddy': u_w_zeddy_teddy,\
       'u_w_zeddy': u_w_zeddy, 'u_w_teddy':u_w_teddy, 'u_w_zeddy2': u_w_zeddy2,\
       'u_w_zmean_teddy': u_w_zmean_teddy}


w_dudp = {'w_dudp_zmean_tmean': w_dudp_zmean_tmean, 'w_dudp_zmean_teddy': w_dudp_zmean_teddy,\
          'w_dudp_zmean': w_dudp_zmean, 'w_dudp_zmean2': w_dudp_zmean2}


EMF = {'EMF_zeddy': EMF_zeddy, 'EMF_zeddy2': EMF_zeddy2, 'EMF_teddy':EMF_teddy,\
       'EMF_zeddy_tmean': EMF_zeddy_tmean, 'EMF_zeddy_teddy': EMF_zeddy_teddy,\
       'EMF_zmean_teddy':EMF_zmean_teddy }

['div_EMF_zmean_teddy',
 'div_EMF_zeddy_teddy',
 'div_EMF_teddy',
 'div_EMF_zeddy_tmean',
 'div_EMF_zeddy',
 'div_EMF_zeddy2']



divEMF = {'div_EMF_zeddy': div_EMF_zeddy, 'div_EMF_zeddy2': div_EMF_zeddy2,\
          'div_EMF_teddy':div_EMF_teddy, 'div_EMF_zeddy_tmean':div_EMF_zeddy_tmean, 'div_EMF_zeddy_teddy':div_EMF_zeddy_teddy,\
          'div_EMF_zmean_teddy': div_EMF_zmean_teddy}


zeta_v = {'zeta_v_zmean_tmean': zeta_v_zmean_tmean, 'zeta_v_zmean_teddy': zeta_v_zmean_teddy,\
          'zeta_v_zmean': zeta_v_zmean, 'zeta_v_zmean2': zeta_v_zmean2, 'zeta_mean': zeta_mean}


du_w_dp = {'du_w_zeddy_tmean_dp': du_w_zeddy_tmean_dp, 'du_w_zeddy_teddy_dp': du_w_zeddy_teddy_dp,\
           'du_w_zeddy_dp': du_w_zeddy_dp, 'du_w_teddy_dp':du_w_teddy_dp, 'du_w_zeddy2_dp': du_w_zeddy2_dp,\
           'du_w_zmean_teddy_dp':du_w_zmean_teddy_dp}


momentum_terms = {'fv_zmean_tmean':fv_zmean_tmean,'du_dt':du_dt, 'w_dudp':w_dudp, 'zeta_v':zeta_v,\
                  'EMF':EMF, 'divEMF':divEMF, 'u_w':u_w, 'du_w_dp':du_w_dp}

coord          = {'lat':latn, 'pres': pres, 'months_per_year': MONTHS}

raw_data1       = {'MSE': MSE, 'CpT':CpT, 'Lq':Lq, 'gZ':gZ, 'mv':mv, 'v':v}

save(destination+"momentum_terms.hkl", momentum_terms)
save(destination+"coord.hkl", coord)
save(destination+"raw_data.hkl", raw_data1)

logging.debug("Saved averaged momentum terms")

end = ti.time()
logging.debug("Looks great !! Time taken --> "+str(end-start))


