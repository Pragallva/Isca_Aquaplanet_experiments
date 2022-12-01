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

edge=str(dirc[1])
land=str(dirc[2])
ocean=str(dirc[3])
num=str(dirc[4]) ## represents an annual year of data

log_directory='/project2/tas1/pragallva/Fall_quarter_2018/codes/shell_script/reverse/log/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'
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

destination='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/reverse/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m/'
make_sure_path_exists(destination);

del_m_mmc=[]
del_m_mmc_sensible=[]
del_m_mmc_pot=[]
del_m_mmc_moist=[]

del_v_mmc=[]
del_v_mmc_sensible=[]
del_v_mmc_pot=[]
del_v_mmc_moist=[]

del_m_mmc_vert=[]
del_m_mmc_sensible_vert=[]
del_m_mmc_pot_vert=[]
del_m_mmc_moist_vert=[]

del_v_mmc_vert=[]
del_v_mmc_sensible_vert=[]
del_v_mmc_pot_vert=[]
del_v_mmc_moist_vert=[]

non_vertical_fluxes=[]
vertical_fluxes=[]

for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/reverse/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+str(i)+'/'
 non_vertical_fluxes.append(load(source+"MMC_thermo_dynamic_dic.hkl"))
 vertical_fluxes.append(load(source+"MMC_thermo_dynamic_dic_vert.hkl"))

 del_m_mmc.append(non_vertical_fluxes[i]['del_m_mmc'])
 del_v_mmc.append(non_vertical_fluxes[i]['del_v_mmc'])
 
 del_m_mmc_moist.append(non_vertical_fluxes[i]['del_m_mmc_moist'])
 del_v_mmc_moist.append(non_vertical_fluxes[i]['del_v_mmc_moist'])

 del_m_mmc_sensible.append(non_vertical_fluxes[i]['del_m_mmc_sensible'])
 del_v_mmc_sensible.append(non_vertical_fluxes[i]['del_v_mmc_sensible'])
 
 del_m_mmc_pot.append(non_vertical_fluxes[i]['del_m_mmc_pot'])
 del_v_mmc_pot.append(non_vertical_fluxes[i]['del_v_mmc_pot'])
  
 del_m_mmc_vert.append(vertical_fluxes[i]['del_m_mmc'])
 del_v_mmc_vert.append(vertical_fluxes[i]['del_v_mmc'])
 
 del_m_mmc_moist_vert.append(vertical_fluxes[i]['del_m_mmc_moist'])
 del_v_mmc_moist_vert.append(vertical_fluxes[i]['del_v_mmc_moist'])
 
 del_m_mmc_sensible_vert.append(vertical_fluxes[i]['del_m_mmc_sensible'])
 del_v_mmc_sensible_vert.append(vertical_fluxes[i]['del_v_mmc_sensible'])

 del_m_mmc_pot_vert.append(vertical_fluxes[i]['del_m_mmc_pot'])
 del_v_mmc_pot_vert.append(vertical_fluxes[i]['del_v_mmc_pot'])

 logging.debug("----"+str(i)+"----")
 
print np.shape( np.mean(del_m_mmc_sensible,axis=0) )
print np.shape( np.mean(del_v_mmc_sensible, axis=0) )

coord        =(load(source+"coord_dic.hkl"))
lat          =coord['lat']
lon          =coord['lon']
no_of_plevels=coord["no_of_plevels"]

def M(x):
    x=np.array(x)
    return np.mean(x,axis=0)

del_m_mmc=M(del_m_mmc); del_v_mmc=M(del_v_mmc);
del_m_mmc_sensible=M(del_m_mmc_sensible); del_v_mmc_sensible=M(del_v_mmc_sensible);
del_m_mmc_pot=M(del_m_mmc_pot); del_v_mmc_pot=M(del_v_mmc_pot);
del_m_mmc_moist=M(del_m_mmc_moist); del_v_mmc_moist=M(del_v_mmc_moist);

del_m_mmc_vert=M(del_m_mmc_vert); del_v_mmc_vert=M(del_v_mmc_vert);
del_m_mmc_sensible_vert=M(del_m_mmc_sensible_vert); del_v_mmc_sensible_vert=M(del_v_mmc_sensible_vert);
del_m_mmc_pot_vert=M(del_m_mmc_pot_vert); del_v_mmc_pot_vert=M(del_v_mmc_pot_vert);
del_m_mmc_moist_vert=M(del_m_mmc_moist_vert); del_v_mmc_moist_vert=M(del_v_mmc_moist_vert);

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS=12
MM_interp       =np.zeros((LATN,MONTHS))
del_m_mmc_interp=np.zeros((LATN,MONTHS))
del_v_mmc_interp=np.zeros((LATN,MONTHS))

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

def zon_int_vert(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None,None]))*a
    return y/10**15

#lat=lat.astype(TE_interp.dtype)
import scipy.integrate as integrate
def integrated(x):
    l=np.deg2rad(latn[:,None])
    
    def A(X):
        Y=(X-np.average(X, axis=0, weights=np.cos(np.deg2rad(latn)))[None,:])
        return Y
    
    x=A(x)
    x=x*np.cos(l)
    int_x  =integrate.cumtrapz(x[::-1,:],l[::-1],axis=0,initial=None) #  (This is basically integration from - 90 deg)
    int_x_r=integrate.cumtrapz(x        ,l      ,axis=0,initial=None) #  (This is basically integration from + 90 deg) 
    avg_int_r=2*np.pi*a**2*(int_x_r[:-1,:])#int_x[::-1,:][1:,:]+/2.0
    return avg_int_r/10**15
              
def interpolate_non_vertical(X):
    interp=np.zeros((LATN,MONTHS))
    for m in range(12):
            interpolation_function = interp1d(lat, X[m,:],kind='linear')
            interp[:,m]=interpolation_function(latn)
    return interp
       
def interpolate_vertical(X):
    interp=np.zeros((LATN,no_of_plevels,MONTHS))
    for m in range(12):
        for p in range(no_of_plevels):
            interpolation_function = interp1d(lat, X[m,p,:],kind='linear')
            interp[:,p,m]=interpolation_function(latn)
    return interp
     
    
del_m_mmc=interpolate_non_vertical(del_m_mmc); 
del_v_mmc=interpolate_non_vertical(del_v_mmc);

del_m_mmc_sensible=interpolate_non_vertical(del_m_mmc_sensible); 
del_v_mmc_sensible=interpolate_non_vertical(del_v_mmc_sensible);

del_m_mmc_pot=interpolate_non_vertical(del_m_mmc_pot); 
del_v_mmc_pot=interpolate_non_vertical(del_v_mmc_pot);

del_m_mmc_moist=interpolate_non_vertical(del_m_mmc_moist); 
del_v_mmc_moist=interpolate_non_vertical(del_v_mmc_moist);



del_m_mmc_vert=interpolate_vertical(del_m_mmc_vert); 
del_v_mmc_vert=interpolate_vertical(del_v_mmc_vert);

del_m_mmc_sensible_vert=interpolate_vertical(del_m_mmc_sensible_vert); del_v_mmc_sensible_vert=interpolate_vertical(del_v_mmc_sensible_vert);

del_m_mmc_pot_vert=interpolate_vertical(del_m_mmc_pot_vert); 
del_v_mmc_pot_vert=interpolate_vertical(del_v_mmc_pot_vert);

del_m_mmc_moist_vert=interpolate_vertical(del_m_mmc_moist_vert); 
del_v_mmc_moist_vert=interpolate_vertical(del_v_mmc_moist_vert);


logging.debug("Calculated interpolation")

############################################
######### Calculate the derivatives ########
############################################

dtheta=np.radians(latn[1]-latn[0])

def spher_div(x):
       N=50
       div=np.copy(x)
       for m in range(12):
           div[:,m]= smooth( np.gradient((x[:,m])*(np.cos(np.radians(latn))),dtheta)/( R*np.cos(np.radians(latn[:])) ),N)
       return div

mmc_decompose_interp = {'del_m_mmc':zon_int(del_m_mmc),'del_v_mmc':zon_int(del_v_mmc),'del_m_mmc_moist':zon_int(del_m_mmc_moist),'del_m_mmc_sensible':zon_int(del_m_mmc_sensible),'del_m_mmc_pot':zon_int(del_m_mmc_pot),'del_v_mmc_moist':zon_int(del_v_mmc_moist),'del_v_mmc_sensible':zon_int(del_v_mmc_sensible),'del_v_mmc_pot':zon_int(del_v_mmc_pot),'latn':latn, } #'KE':zon_int(KE_interp) 


save(destination+'mmc_decompose_interp.hkl', mmc_decompose_interp)
logging.debug("Saved decompose mmc dictionary")

mmc_decompose_vert_interp = {'del_m_mmc':zon_int_vert(del_m_mmc_vert),'del_v_mmc':zon_int_vert(del_v_mmc_vert),'del_m_mmc_moist':zon_int_vert(del_m_mmc_moist_vert),'del_m_mmc_sensible':zon_int_vert(del_m_mmc_sensible_vert),'del_m_mmc_pot':zon_int_vert(del_m_mmc_pot_vert),'del_v_mmc_moist':zon_int_vert(del_v_mmc_moist_vert),'del_v_mmc_sensible':zon_int_vert(del_v_mmc_sensible_vert),'del_v_mmc_pot':zon_int_vert(del_v_mmc_pot_vert),'latn':latn } #'KE':zon_int(KE_interp) 

save(destination+'mmc_decompose_vert_interp.hkl', mmc_decompose_vert_interp)
logging.debug("Saved decompose mmc dictionary")


logging.debug("Looks great !!")


