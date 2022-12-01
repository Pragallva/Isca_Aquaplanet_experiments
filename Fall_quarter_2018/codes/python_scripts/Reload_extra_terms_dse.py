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

log_directory='/project2/tas1/pragallva/Fall_quarter_2018/codes/shell_script/log/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'
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

destination='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m/'
make_sure_path_exists(destination);

tendency=[]

dhdt=[]
dsensdt=[]
dpotdt=[]
dmoistdt=[]

dhdt_v=[]
dsensdt_v=[]
dpotdt_v=[]
dmoistdt_v=[]

for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+str(i)+'/'
 tendency.append(load(source+"tendency_Wm2_dic.hkl"))
    
 dhdt.append(tendency[i]['dhdt'])
 dsensdt.append(tendency[i]['dsensdt'])
 dmoistdt.append(tendency[i]['dmoistdt'])
 dpotdt.append(tendency[i]['dpotdt'])

 dhdt_v.append(tendency[i]['dhdt_vert'])
 dsensdt_v.append(tendency[i]['dsensdt_vert'])
 dmoistdt_v.append(tendency[i]['dmoistdt_vert'])
 dpotdt_v.append(tendency[i]['dpotdt_vert'])

 logging.debug("----"+str(i)+"----")
 
coord=(load(source+"coord_dic.hkl"))
lat      =coord['lat']
lon      =coord['lon']
no_of_plevels=coord["no_of_plevels"]

def M(x):
    return np.mean(np.array(x),axis=0)

dhdt    = M(dhdt);     #dhdt_v     = M(dhdt_v); 
dsensdt = M(dsensdt);  dsensdt_v  = M(dsensdt_v); 
dmoistdt= M(dmoistdt); dmoistdt_v = M(dmoistdt_v); 
dpotdt  = M(dpotdt);   dpotdt_v   = M(dpotdt_v); 

LAT=len(lat)
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS=12

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

lat=lat.astype(np.float64)
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

def integratedv(x):
    l=np.deg2rad(latn[:,None,None])

    def A(X):
        Y=(X-np.average(X, axis=0, weights=np.cos(np.deg2rad(latn)))[None,...])
        return Y

    x=A(x)
    x=x*np.cos(l)
    int_x  =integrate.cumtrapz(x[::-1,:,:],l[::-1],axis=0,initial=None) #  (This is basically integration from - 90 deg)
    int_x_r=integrate.cumtrapz(x        ,l      ,axis=0,initial=None) #  (This is basically integration from + 90 deg) 
    avg_int_r=2*np.pi*a**2*(int_x_r[:-1,:,:])#int_x[::-1,:][1:,:]+/2.0
    return avg_int_r/10**15
 
def interpolatev(X):
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

dhdt    = interpolate(dhdt);     #dhdt_v     = interpolatev(dhdt_v);
dsensdt = interpolate(dsensdt);  dsensdt_v  = interpolatev(dsensdt_v);
dmoistdt= interpolate(dmoistdt); dmoistdt_v = interpolatev(dmoistdt_v);
dpotdt  = interpolate(dpotdt);   dpotdt_v   = interpolatev(dpotdt_v);

dhdt_v     = dsensdt_v+dmoistdt_v

logging.debug("Calculated interpolation")

############################################
######### Calculate the derivatives ########
############################################

dtheta=np.radians(latn[1]-latn[0])

tendency_Wm2_zonal_dic = {"dhdt":dhdt,"dhdt_vert":dhdt_v,'dsensdt':dsensdt,'dsensdt_vert':dsensdt_v, 'dmoistdt':dmoistdt,'dmoistdt_vert': dmoistdt_v, 'dpotdt':dpotdt, 'dpotdt_vert':dpotdt_v}

tendency_PW_zonal_dic = {"dhdt":integrated(dhdt),"dhdt_vert":integratedv(dhdt_v),'dsensdt':integrated(dsensdt),'dsensdt_vert':integratedv(dsensdt_v), 'dmoistdt':integrated(dmoistdt),'dmoistdt_vert': integratedv(dmoistdt_v), 'dpotdt':integrated(dpotdt), 'dpotdt_vert':integratedv(dpotdt_v)}

save(destination+'xtendency_Wm2_zonal_dic.hkl', tendency_Wm2_zonal_dic)
save(destination+'xtendency_PW_zonal_dic.hkl', tendency_PW_zonal_dic)
logging.debug("Saved tendency interp dictionary")

logging.debug("Looks great !!")


