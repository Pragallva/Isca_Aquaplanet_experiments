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
log_directory='/project2/tas1/pragallva/Fall_quarter_2017/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

# try:
#    import cPickle as pickle        Unfortunately pickle doesn't work
# except:
#    import pickle

logging.debug("** INTERPOLATION **")

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

            
source='/project2/tas1/pragallva/Fall_quarter_2017/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'

decomposed_fluxes=load(source+"zonal_decomposed_fluxes_dic.hkl")
coord=load(source+"coord_dic.hkl")
radiation=load(source+"zonal_radiation_dic.hkl")

TE       =decomposed_fluxes['TE_flux']
SE       =decomposed_fluxes['SE_flux']
MM       =decomposed_fluxes['MM_flux']
KE       =decomposed_fluxes['KE_flux']

toa      =radiation['TOA_d']
sfc      =radiation['SFC_u']
dhdt     =decomposed_fluxes['dhdt']
SWABS    =radiation['SWABS']
SHF      =radiation['SHF']
OLR      =radiation['olr']

lat      =coord['lat']
lon      =coord['lon']

LAT=len(lat)
MSE=MM+SE+TE+KE
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1)
LATN=len(latn)   

MONTHS=12
TE_interp=np.zeros((LATN,MONTHS))
SE_interp=np.zeros((LATN,MONTHS))
MM_interp=np.zeros((LATN,MONTHS))
KE_interp=np.zeros((LATN,MONTHS))
MSE_interp=np.zeros((LATN,MONTHS))
toa_interp=np.zeros((LATN,MONTHS))
sfc_interp=np.zeros((LATN,MONTHS))
dhdt_interp=np.zeros((LATN,MONTHS))
SWABS_interp=np.zeros((LATN,MONTHS))
SHF_interp=np.zeros((LATN,MONTHS))
OLR_interp=np.zeros((LATN,MONTHS))

def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

import scipy.integrate as integrate
def integrated(x):
    l=np.deg2rad(latn[:,None])
    x=x*np.cos(l)
    int_x  =integrate.cumtrapz(x[::-1,:],l[::-1],axis=0,initial=None) #  (This is basically integration from - 90 deg)
    int_x_r=integrate.cumtrapz(x        ,l      ,axis=0,initial=None) #  (This is basically integration from + 90 deg) 
    avg_int_r=2*np.pi*a**2*(int_x[::-1,:][1:,:]+int_x_r[:-1,:])/2.0
    return avg_int_r/10**15

for m in range(12):

        MM_interpolation_function = interp1d(lat, MM[m,:],kind='cubic')
        MM_interp[:,m]=MM_interpolation_function(latn)

        TE_interpolation_function = interp1d(lat, TE[m,:],kind='cubic')
        TE_interp[:,m]=TE_interpolation_function(latn)

        SE_interpolation_function= interp1d(lat, SE[m,:],kind='cubic')
        SE_interp[:,m]= SE_interpolation_function(latn)
        
        MSE_interpolation_function= interp1d(lat, MSE[m,:],kind='cubic')
        MSE_interp[:,m]= MSE_interpolation_function(latn)
        
        KE_interpolation_function= interp1d(lat, KE[m,:],kind='cubic')
        KE_interp[:,m]= KE_interpolation_function(latn)
        
        toa_interpolation_function= interp1d(lat, toa[m,:],kind='cubic')
        toa_interp[:,m]= toa_interpolation_function(latn)
        
        sfc_interpolation_function= interp1d(lat, sfc[m,:],kind='cubic')
        sfc_interp[:,m]= sfc_interpolation_function(latn)
        
        dhdt_interpolation_function= interp1d(lat, dhdt[m,:],kind='cubic')
        dhdt_interp[:,m]= dhdt_interpolation_function(latn)
        
        SWABS_interpolation_function= interp1d(lat, SWABS[m,:],kind='cubic')
        SWABS_interp[:,m]= SWABS_interpolation_function(latn)
        
        SHF_interpolation_function= interp1d(lat, SHF[m,:],kind='cubic')
        SHF_interp[:,m]= SHF_interpolation_function(latn)
        
        OLR_interpolation_function= interp1d(lat, OLR[m,:],kind='cubic')
        OLR_interp[:,m]= OLR_interpolation_function(latn)
        
logging.debug("Calculated interpolation")
        
flux_interp_dict={'MM':zon_int(MM_interp),'SE':zon_int(SE_interp),'TE':zon_int(TE_interp),'MSE':zon_int(MSE_interp),'KE':zon_int(KE_interp),'TOA_d':integrated(toa_interp),'SFC_u':integrated(sfc_interp),'dhdt':integrated(dhdt_interp),'SWABS':integrated(SWABS_interp),'SHF':integrated(SHF_interp),'olr':integrated(OLR_interp),'latn':latn,'latnr':latn[1:-1]}        

############################################
######### Calculate the derivatives ########
############################################

dtheta=np.radians(latn[1]-latn[0])

def spher_div(x):
       N=10
       div=np.copy(x)
       for m in range(12):
           div[:,m]= smooth( np.gradient((x[:,m])*(np.cos(np.radians(latn))),dtheta)/( R*np.cos(np.radians(latn[:])) ),N)
       return div
           
div_TE_flux=spher_div(TE_interp)
div_SE_flux=spher_div(SE_interp)
div_MM_flux=spher_div(MM_interp)
div_KE_flux=spher_div(KE_interp)
div_MSE_flux=spher_div(MSE_interp)
div_NE_flux=toa_interp+sfc_interp-dhdt_interp

logging.debug("Calculated divergence")

div_flux_dict={'MM':div_MM_flux,'SE':div_SE_flux,'TE':div_TE_flux,'MSE':div_MSE_flux,'KE':div_KE_flux,'TOA_d':toa_interp,'SFC_u':sfc_interp,'dhdt':dhdt_interp,'SWABS':SWABS_interp,'SHF':SHF_interp,'olr':OLR_interp,'latn':latn,'latnr':latn[1:-1]} 

save(source+'flux_interp_dict.hkl', flux_interp_dict)
logging.debug("Saved flux interp dictionary")
save(source+'div_flux_dict.hkl', div_flux_dict)
logging.debug("Saved flux divergence interp")

logging.debug("Looks great !!")


