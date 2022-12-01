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

log_directory='/project2/tas1/pragallva/Fall_quarter_2018/codes/shell_script/log/'+'baroclinic'+edge+'_la'+land+'m_oc'+ocean+'m'
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

source1='/project2/tas1/pragallva/Spring_quarter_2018/exp_data/aqua_isca100m/
one_year="aqua_isca100m0.nc"
ncfile=source1+one_year
lat=nc.Dataset(ncfile,'r').variables['lat'][:]

stream_data=[]
uv_data=[]
Z_T_Q_data=[]
T=[]
Z=[]
U=[]


for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/data_in_pres_coord/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+str(i)+'/'
 Z_T_Q_data.append(load(source+"ht_temp_sphum.hkl"))   
 uv_data.append(load(source+"u_v.hkl"))   
 T.append(Z_T_Q_data[i]['temp'])
 Z.append(Z_T_Q_data[i]['ht'])
 U.append(uv_data[i]['u'])

 logging.debug("----"+str(i)+"----T, U, Z ")
     
def M(x):
    return np.mean(x,axis=0)

T=M(T) ; Z=M(Z); U=M(U)

Pres = np.array([0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0])
no_of_plevels = len(Pres)
P_half=np.append(0,Pres)

logging.debug(" Loaded T U Z P")

Rk=287 # --> Gas constant
Cp= 1004.64 # J/kg/deg
g = 9.8
L = 2.500e6   # J/kg
omega = 7.2921*10**-5 # ---> Rotation rate of earth rad/s
#f = 2*omega*np.sin(np.deg2rad(lat_a)) # --> Coriolis frequency
a=6371.0e3
R=a

logging.debug("Defined constants")

#### Interpolation ####

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  
        
def interpolate_plevels(X):
    interp=np.zeros((LATN,no_of_plevels,X.shape[0]))
    for m in range(X.shape[0]):
        for p in range(no_of_plevels):
            interpolation_function = interp1d(lat, X[m,p,:],kind='linear')
            interp[:,p,m]=interpolation_function(latn)
    return interp


T_interp =interpolate_plevels(T)  #lat, plev, time
logging.debug("Interpolate T")
U_interp =interpolate_plevels(U)  #lat, plev, time
logging.debug("Interpolate U")
Z_interp =interpolate_plevels(Z)  #lat, plev, time
logging.debug("Interpolate Z")

theta = T_interp*((Pres/1000.0)**(-Rk/Cp))[None,:,None] ## Potential temperature

logging.debug("Calculated theta")

def P(p):
    mini=np.min(np.abs(Pres-p))
    y=np.squeeze(np.where( np.abs(Pres-p)==mini  ) )
    return y

def weighted(arg):
    weights= (P_half[1:]-P_half[:-1])[None,:,None] 
    w      = arg*weights
    return w

pmin = 700.0; pmax=850.0
pmean       = lambda pmin, pmax: (pmin+pmax)/2.0
deltap      = lambda pmin, pmax: (pmax-pmin)
p00         = 1000.0
                
delta_theta = lambda pmin, pmax: theta[:,P(pmax),:]-theta[:,P(pmin),:]
           
def theta_avg(pmin, pmax):
    y          = weighted(theta)
    integrated = y[:,P(pmin):P(pmax),:].sum(axis=1)
    delta_p    = Pres[P(pmax)]-Pres[P(pmin)]
    return integrated/delta_p # (lat, time)

def spher_div(x):
       dtheta=np.deg2rad(latn[1]-latn[0])
       N=50
       div=np.copy(x)
       for t in range( (x.shape[-1]) ):
          div[:,t]= smooth( np.gradient((x[:,t]),dtheta)/( R*np.cos(np.radians(latn)) ),N) 
       return div # (lat, time)
           
def vert_grad(x):
       N=1
       div=np.copy(x)
       for t in range( (x.shape[-1]) ):
           for l in range( (x.shape[0]) ):
              div[l,:,t]= smooth( np.gradient((x[l,:,t]),Pres),N)
       return div


N_sq_sfc = -(g**2/Rk)* (pmean(pmin,pmax)/deltap(pmin,pmax)) * (p00/pmean(pmin,pmax))**(Rk/Cp) * (delta_theta(pmin,pmax))/(theta_avg(pmin,pmax)) # (lat, time)
logging.debug("Calculated N_sq_sfc")


dtheta_dy    = spher_div( theta_avg(pmin,pmax) ) # (lat, time)
logging.debug("Calculated dtheta_dy_sfc")

fo_dUdz_sfc  = -1*g*dtheta_dy/theta_avg(pmin,pmax)  ## Have not divided by f0 (lat, time)
logging.debug("Calculated fo_dUdz_sfc")
           
sfc_baroclinicity = {'N_sq':N_sq_sfc,'fo_dUdz': fo_dUdz_sfc, 'latn': latn } # (lat, time)
logging.debug("Dictionary surface baroclinicity")

### From formula  ####
dUdp      = vert_grad(U_interp)  # (lat, plev, time)
dZdp      = vert_grad(Z_interp)  # (lat, plev, time)
dUdZ      = dUdp/dZdp            # (lat, plev, time)
logging.debug("Calculated dUdz")
           
dthetadp = vert_grad(theta)  
N_sq     = -(g**2*Pres[None,:,None]/(R*T_interp))* dthetadp/(theta)
logging.debug("Calculated N_sq")
                  
baroclinicity = {'N_sq':N_sq,'dUdz': dUdZ, 'Pres': Pres, 'latn': latn }       
logging.debug("Dictionary baroclinicity")

destination='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/data_in_pres_coord/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m/'
make_sure_path_exists(destination)

save(destination+"baroclinicity_sfc.hkl", sfc_baroclinicity)
save(destination+"baroclinicity.hkl", baroclinicity)
           
logging.debug("Saved baroclinicity")
logging.debug("Looks great !!")


