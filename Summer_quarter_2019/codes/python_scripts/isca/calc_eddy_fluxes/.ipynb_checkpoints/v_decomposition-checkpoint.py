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

logging.basicConfig( filename = log_directory+'v_decomp'+exp_dir+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
         
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


def df(y, deno=0):
    rot_rate  = 7.29*10**-5
    f=(2*rot_rate*np.sin(np.deg2rad(lat)))[None,None,None,:]+np.zeros(np.shape(y))
    if deno==0:
       y1=y/f
    else :
       y1=y/(f+np.array(zeta_mean))
    return y1
    
def M1(x):
    return np.array(x)

v           = df(M1(fv_zmean_tmean),   0)

v_dudt      = df(+M1(du_dt),   0)
v_zetaV     = df(-M1(zeta_v_zmean_tmean),   0)
v_wdudp     = df(+M1(w_dudp_zmean_tmean),  0)
v_divEMF_SE = df(+M1(div_EMF_zeddy_tmean), 0)
v_divEMF_TE = df(+M1(div_EMF_teddy), 0)
v_du_wdp_SE = df(+M1(du_w_zeddy_tmean_dp), 0)
v_du_wdp_TE = df(+M1(du_w_teddy_dp), 0)


v_dudt1      = df(+M1(du_dt),   1)
v_wdudp1     = df(-M1(w_dudp_zmean_tmean),  1)
v_divEMF_SE1 = df(+M1(div_EMF_zeddy_tmean), 1)
v_divEMF_TE1 = df(+M1(div_EMF_teddy), 1)
v_du_wdp_SE1 = df(+M1(du_w_zeddy_tmean_dp), 1)
v_du_wdp_TE1 = df(+M1(du_w_teddy_dp), 1)


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
         


v_decomp  = {'v_dudt': interpolate_vertical(M(v_dudt)), 'v_zetaV': interpolate_vertical(M(v_zetaV)),\
             'v_wdudp': interpolate_vertical(M(v_wdudp)), 'v_divEMF_SE': interpolate_vertical(M(v_divEMF_SE)),\
             'v_divEMF_TE': interpolate_vertical(M(v_divEMF_TE)), 'v_du_wdp_SE':  interpolate_vertical(M(v_du_wdp_SE)),\
             'v_du_wdp_TE': interpolate_vertical(M(v_du_wdp_TE)), 'v': interpolate_vertical(M(v))}

v_decomp1  = {'v_dudt1': interpolate_vertical(M(v_dudt1)),'v1':interpolate_vertical(M(v)),\
             'v_wdudp1': interpolate_vertical(M(v_wdudp1)), 'v_divEMF_SE1': interpolate_vertical(M(v_divEMF_SE1)),\
             'v_divEMF_TE1': interpolate_vertical(M(v_divEMF_TE1)), 'v_du_wdp_SE1':  interpolate_vertical(M(v_du_wdp_SE1)),\
             'v_du_wdp_TE1': interpolate_vertical(M(v_du_wdp_TE1))}


save(destination+"v_decomp.hkl", {'v_decomp':v_decomp,'v_decomp1':v_decomp1})
logging.debug("Saved averaged momentum terms")

end = ti.time()
logging.debug("Looks great !! Time taken --> "+str(end-start))


