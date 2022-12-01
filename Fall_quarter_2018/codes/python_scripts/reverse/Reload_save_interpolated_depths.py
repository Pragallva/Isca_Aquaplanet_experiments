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

decomposed_fluxes=[]
decomposed_fluxes_vert=[]
radiation=[]
raw_data=[]
raw_MSE=[]

TE=[]
SE=[]
MM=[]
#KE=[]

TE_sensible=[]
TE_moist=[]
TE_pot=[]

SE_sensible=[]
SE_moist=[]
SE_pot=[]

MM_sensible=[]
MM_moist=[]
MM_pot=[]

TE_vert=[]
SE_vert=[]
MM_vert=[]

toa=[]
sfc=[]
dhdt=[]
SWABS=[]
SHF=[]
OLR=[]

SW_toa_d=[]
shflx_u=[]
lhflx_u=[]
SW_sfc_d=[]
LW_sfc_d=[]
precip=[]

T=[]
U=[]
V=[]
Z=[]
q=[]

all_MSE=[]


for i in range(0,int(num)):
 source='/project2/tas1/pragallva/Fall_quarter_2018/post_process_data/reverse/'+'HC'+edge+'_la'+land+'m_oc'+ocean+'m'+str(i)+'/'
 decomposed_fluxes.append(load(source+"zonal_decomposed_fluxes_dic.hkl"))
 decomposed_fluxes_vert.append(load(source+"zonal_decomposed_fluxes_dic_vert.hkl"))
 radiation.append(load(source+"zonal_radiation_dic.hkl"))
 raw_data.append(load(source+"T_uv_Z.hkl"))
 raw_MSE.append(load(source+"all_MSE_flux.hkl"))
    
 T.append(raw_data[i]['temp'])
 U.append(raw_data[i]['u'])
 V.append(raw_data[i]['v'])
 Z.append(raw_data[i]['Z'])
 q.append(raw_data[i]['q'])
 all_MSE.append(raw_MSE[i]['MSE_flux'])

 TE.append(decomposed_fluxes[i]['TE_flux'])
 SE.append(decomposed_fluxes[i]['SE_flux'])
 MM.append(decomposed_fluxes[i]['MM_flux'])
 #KE.append(decomposed_fluxes[i]['KE_flux'])

 TE_sensible.append(decomposed_fluxes[i]['TE_sensible'])
 SE_sensible.append(decomposed_fluxes[i]['SE_sensible'])
 MM_sensible.append(decomposed_fluxes[i]['MM_sensible'])
    
 TE_moist.append(decomposed_fluxes[i]['TE_moist'])
 SE_moist.append(decomposed_fluxes[i]['SE_moist'])
 MM_moist.append(decomposed_fluxes[i]['MM_moist'])

 #TE_pot.append(decomposed_fluxes[i]['TE_pot'])
 SE_pot.append(decomposed_fluxes[i]['SE_pot'])
 MM_pot.append(decomposed_fluxes[i]['MM_pot'])
    
 TE_vert.append(decomposed_fluxes_vert[i]['TE_flux'])
 SE_vert.append(decomposed_fluxes_vert[i]['SE_flux'])
 MM_vert.append(decomposed_fluxes_vert[i]['MM_flux'])

 toa.append(radiation[i]['TOA_d'])
 sfc.append(radiation[i]['SFC_u'])
 dhdt.append(decomposed_fluxes[i]['dhdt'])
 SWABS.append(radiation[i]['SWABS'])
 SHF.append(radiation[i]['SHF'])
 OLR.append(radiation[i]['olr'])
    
 SW_toa_d.append(radiation[i]['SW_toa_d'])
 SW_sfc_d.append(radiation[i]['SW_sfc_d'])
 LW_sfc_d.append(radiation[i]['LW_sfc_d'])
 shflx_u.append(radiation[i]['shflx_u'])
 lhflx_u.append(radiation[i]['lhflx_u'])
 precip.append(radiation[i]['precip'])   
 logging.debug("----"+str(i)+"----")
 
coord=(load(source+"coord_dic.hkl"))
lat      =coord['lat']
lon      =coord['lon']
no_of_plevels=coord["no_of_plevels"]

def M(x):
    return np.mean(x,axis=0)

TE=M(TE) ; SE=M(SE); MM=M(MM); #KE=M(KE); 
TE_sensible=M(TE_sensible) ; TE_moist=M(TE_moist); 
SE_sensible=M(SE_sensible) ; SE_moist=M(SE_moist); SE_pot=M(SE_pot);
MM_sensible=M(MM_sensible) ; MM_moist=M(MM_moist); MM_pot=M(MM_pot)
TE_pot = TE - (TE_sensible + TE_sensible + TE_moist)


TE_vert=M(TE_vert); SE_vert=M(SE_vert); MM_vert=M(MM_vert)

toa=M(toa); sfc=M(sfc); dhdt=M(dhdt); SWABS=M(SWABS); SHF=M(SHF); OLR=M(OLR)
T=M(T); U=M(U); V=M(V); Z=M(Z); q=M(q); all_MSE=M(all_MSE)
SW_toa_d=M(SW_toa_d) ; SW_sfc_d=M(SW_sfc_d); LW_sfc_d=M(LW_sfc_d); shflx_u=M(shflx_u); lhflx_u=M(lhflx_u)
precip=M(precip)

LAT=len(lat)
MSE=MM+SE+TE#+KE
a=6371.0e3
R=a

latn=np.arange(-87.0,87.1,0.1).astype(np.float64)
LATN=len(latn)  

MONTHS=12
TE_interp=np.zeros((LATN,MONTHS))
SE_interp=np.zeros((LATN,MONTHS))
MM_interp=np.zeros((LATN,MONTHS))
#KE_interp=np.zeros((LATN,MONTHS),dtype=float)
MSE_interp=np.zeros((LATN,MONTHS))
toa_interp=np.zeros((LATN,MONTHS))
sfc_interp=np.zeros((LATN,MONTHS))
dhdt_interp=np.zeros((LATN,MONTHS))
SWABS_interp=np.zeros((LATN,MONTHS))
SHF_interp=np.zeros((LATN,MONTHS))
OLR_interp=np.zeros((LATN,MONTHS))
precip_interp=np.zeros((LATN,MONTHS))

T_interp=np.zeros((LATN,no_of_plevels,MONTHS))
U_interp=np.copy(T_interp)
V_interp=np.copy(T_interp)
q_interp=np.copy(T_interp)
Z_interp=np.copy(T_interp)
all_MSE_interp=np.copy(T_interp)

TE_vert_interp=np.copy(T_interp)
SE_vert_interp=np.copy(T_interp)
MM_vert_interp=np.copy(T_interp)

SW_toa_d_interp=np.copy(TE_interp)
SW_sfc_d_interp=np.copy(TE_interp)
LW_sfc_d_interp=np.copy(TE_interp)
shflx_u_interp=np.copy(TE_interp)
lhflx_u_interp=np.copy(TE_interp)

TE_sens_interp=np.zeros((LATN,MONTHS))
SE_sens_interp=np.zeros((LATN,MONTHS))
MM_sens_interp=np.zeros((LATN,MONTHS))

TE_moist_interp=np.zeros((LATN,MONTHS))
SE_moist_interp=np.zeros((LATN,MONTHS))
MM_moist_interp=np.zeros((LATN,MONTHS))

TE_pot_interp=np.zeros((LATN,MONTHS))
SE_pot_interp=np.zeros((LATN,MONTHS))
MM_pot_interp=np.zeros((LATN,MONTHS))


def zon_int(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None]))*a
    return y/10**15

lat=lat.astype(TE_interp.dtype)
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


for m in range(12):

        MM_interpolation_function = interp1d(lat, MM[m,:],kind='linear')
        MM_interp[:,m]=MM_interpolation_function(latn)

        TE_interpolation_function = interp1d(lat, TE[m,:],kind='linear')
        TE_interp[:,m]=TE_interpolation_function(latn)

        SE_interpolation_function= interp1d(lat, SE[m,:],kind='linear')
        SE_interp[:,m]= SE_interpolation_function(latn)
 

        
        MM_sens_interpolation_function = interp1d(lat, MM_sensible[m,:],kind='linear')
        MM_sens_interp[:,m]=MM_sens_interpolation_function(latn)

        SE_sens_interpolation_function = interp1d(lat, SE_sensible[m,:],kind='linear')
        SE_sens_interp[:,m]=SE_sens_interpolation_function(latn)

        TE_sens_interpolation_function = interp1d(lat, TE_sensible[m,:],kind='linear')
        TE_sens_interp[:,m]=TE_sens_interpolation_function(latn)
        

        
        MM_moist_interpolation_function = interp1d(lat, MM_moist[m,:],kind='linear')
        MM_moist_interp[:,m]=MM_moist_interpolation_function(latn)

        SE_moist_interpolation_function = interp1d(lat, SE_moist[m,:],kind='linear')
        SE_moist_interp[:,m]=SE_moist_interpolation_function(latn)

        TE_moist_interpolation_function = interp1d(lat, TE_moist[m,:],kind='linear')
        TE_moist_interp[:,m]=TE_moist_interpolation_function(latn)
 


        MM_pot_interpolation_function = interp1d(lat, MM_pot[m,:],kind='linear')
        MM_pot_interp[:,m]=MM_pot_interpolation_function(latn)

        SE_pot_interpolation_function = interp1d(lat, SE_pot[m,:],kind='linear')
        SE_pot_interp[:,m]=SE_pot_interpolation_function(latn)

        TE_pot_interpolation_function = interp1d(lat, TE_pot[m,:],kind='linear')
        TE_pot_interp[:,m]=TE_pot_interpolation_function(latn)
        
  

        MSE_interpolation_function= interp1d(lat, MSE[m,:],kind='linear')
        MSE_interp[:,m]= MSE_interpolation_function(latn)
        
        #KE_interpolation_function= interp1d(lat, KE[m,:],kind='cubic')
        #KE_interp[:,m]= KE_interpolation_function(latn)
        
        toa_interpolation_function= interp1d(lat, toa[m,:],kind='linear')
        toa_interp[:,m]= toa_interpolation_function(latn)
        
        sfc_interpolation_function= interp1d(lat, sfc[m,:],kind='linear')
        sfc_interp[:,m]= sfc_interpolation_function(latn)
        
        dhdt_interpolation_function= interp1d(lat, dhdt[m,:],kind='linear')
        dhdt_interp[:,m]= dhdt_interpolation_function(latn)
        
        SWABS_interpolation_function= interp1d(lat, SWABS[m,:],kind='linear')
        SWABS_interp[:,m]= SWABS_interpolation_function(latn)
        
        SHF_interpolation_function= interp1d(lat, SHF[m,:],kind='linear')
        SHF_interp[:,m]= SHF_interpolation_function(latn)
        
        OLR_interpolation_function= interp1d(lat, OLR[m,:],kind='linear')
        OLR_interp[:,m]= OLR_interpolation_function(latn)
        
        SW_toa_d_interpolation_function= interp1d(lat, SW_toa_d[m,:],kind='linear')
        SW_toa_d_interp[:,m]= SW_toa_d_interpolation_function(latn)
        
        SW_sfc_d_interpolation_function= interp1d(lat, SW_sfc_d[m,:],kind='linear')
        SW_sfc_d_interp[:,m]= SW_sfc_d_interpolation_function(latn)
        
        LW_sfc_d_interpolation_function= interp1d(lat, LW_sfc_d[m,:],kind='linear')
        LW_sfc_d_interp[:,m]= LW_sfc_d_interpolation_function(latn)
        
        shflx_u_interpolation_function= interp1d(lat, shflx_u[m,:],kind='linear')
        shflx_u_interp[:,m]= shflx_u_interpolation_function(latn)
        
        lhflx_u_interpolation_function= interp1d(lat, lhflx_u[m,:],kind='linear')
        lhflx_u_interp[:,m]= lhflx_u_interpolation_function(latn)

        precip_interpolation_function=interp1d(lat, precip[m,:],kind='linear')
        precip_interp[:,m]= precip_interpolation_function(latn)
        
def interpolate_raw_data(X):
    interp=np.zeros((LATN,no_of_plevels,MONTHS))
    for m in range(12):
        for p in range(no_of_plevels):
            interpolation_function = interp1d(lat, X[m,p,:],kind='linear')
            interp[:,p,m]=interpolation_function(latn)
    return interp

T_interp      =interpolate_raw_data(T)
Z_interp      =interpolate_raw_data(Z)
U_interp      =interpolate_raw_data(U)
V_interp      =interpolate_raw_data(V)
q_interp      =interpolate_raw_data(q)
all_MSE_interp=interpolate_raw_data(all_MSE)
TE_vert_interp      =interpolate_raw_data(TE_vert)
SE_vert_interp      =interpolate_raw_data(SE_vert)
MM_vert_interp      =interpolate_raw_data(MM_vert)
NE_vert_interp      =interpolate_raw_data(TE_vert+SE_vert+MM_vert)

       
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
           
div_TE_flux=spher_div(TE_interp)
div_SE_flux=spher_div(SE_interp)
div_MM_flux=spher_div(MM_interp)
#div_KE_flux=spher_div(KE_interp)
div_MSE_flux=spher_div(MSE_interp)
div_NE_flux=toa_interp+sfc_interp-dhdt_interp


div_TE_sensible=spher_div(TE_sens_interp)
div_SE_sensible=spher_div(SE_sens_interp)
div_MM_sensible=spher_div(MM_sens_interp)

div_TE_moist=spher_div(TE_moist_interp)
div_SE_moist=spher_div(SE_moist_interp)
div_MM_moist=spher_div(MM_moist_interp)

div_TE_pot=spher_div(TE_pot_interp)
div_SE_pot=spher_div(SE_pot_interp)
div_MM_pot=spher_div(MM_pot_interp)


logging.debug("Calculated divergence")


def spher_div(x):
       N=50
       div=np.copy(x)
       for p in range(no_of_plevels):
           for m in range(12):
               div[:,p,m]= smooth( np.gradient((x[:,p,m])*(np.cos(np.radians(latn))),dtheta)/( R*np.cos(np.radians(latn[:])) ),N)
       return div
           
div_TE_vert_flux=spher_div(TE_vert_interp)
div_SE_vert_flux=spher_div(SE_vert_interp)
div_MM_vert_flux=spher_div(MM_vert_interp)
#div_KE_flux=spher_div(KE_interp)
div_NE_vert_flux=div_TE_vert_flux+div_SE_vert_flux+div_MM_vert_flux



flux_interp_dict={'MM':zon_int(MM_interp),'SE':zon_int(SE_interp),'TE':zon_int(TE_interp),'MM_sens':zon_int(MM_sens_interp),'SE_sens':zon_int(SE_sens_interp),'TE_sens':zon_int(TE_sens_interp),'MM_moist':zon_int(MM_moist_interp),'SE_moist':zon_int(SE_moist_interp),'TE_moist':zon_int(TE_moist_interp),'MM_pot':zon_int(MM_pot_interp),'SE_pot':zon_int(SE_pot_interp),'TE_pot':zon_int(TE_moist_interp),'MSE':zon_int(MSE_interp),'TOA_d':integrated(toa_interp),'SFC_u':integrated(sfc_interp),'dhdt':integrated(dhdt_interp),'SWABS':integrated(SWABS_interp),'SHF':integrated(SHF_interp),'olr':integrated(OLR_interp),'latn':latn, 'SW_toa_d':integrated(SW_toa_d_interp), 'SW_sfc_d':integrated(SW_sfc_d_interp), 'LW_sfc_d':integrated(LW_sfc_d_interp), 'shflx_u':integrated(shflx_u_interp), 'lhflx_u':integrated(lhflx_u_interp),'precip':precip_interp ,'latnr':latn[1:-1]} #'KE':zon_int(KE_interp) 

div_flux_dict={'MM':div_MM_flux,'SE':div_SE_flux,'TE':div_TE_flux,'MM_sens':div_MM_sensible,'SE_sens':div_SE_sensible,'TE_sens':div_TE_sensible,'MM_moist':div_MM_moist,'SE_moist':div_SE_moist,'TE_moist':div_TE_moist,'MM_pot':div_MM_pot,'SE_pot':div_SE_pot,'TE_pot':div_TE_pot,'MSE':div_MSE_flux,'TOA_d':toa_interp,'SFC_u':sfc_interp,'dhdt':dhdt_interp,'SWABS':SWABS_interp,'SHF':SHF_interp,'olr':OLR_interp,'latn':latn,'latnr':latn[1:-1], 'SW_toa_d':(SW_toa_d_interp), 'SW_sfc_d':(SW_sfc_d_interp), 'LW_sfc_d':(LW_sfc_d_interp), 'shflx_u':(shflx_u_interp), 'lhflx_u':(lhflx_u_interp)} ;#,'KE':div_KE_flux

def zon_int_vert(x):
    y=x*2*np.pi*np.cos(np.deg2rad(latn[:,None,None]))*a
    return y/10**15

flux_vert_interp_dict={'MM':zon_int_vert(MM_vert_interp),'SE':zon_int_vert(SE_vert_interp),'TE':zon_int_vert(TE_vert_interp), 'NE':zon_int_vert(NE_vert_interp)} 

div_vert_interp_dict={'MM':zon_int_vert(div_MM_vert_flux),'SE':zon_int_vert(div_SE_vert_flux),'TE':zon_int_vert(div_TE_vert_flux), 'NE':zon_int_vert(div_NE_vert_flux)} 


#'KE':zon_int(KE_interp) 

raw_data_dict={'T':T_interp,'Z':Z_interp,'U':U_interp,'V':V_interp,'q':q_interp,'MSE':all_MSE_interp, 'latn':latn}
save(destination+'raw_data_dict.hkl', raw_data_dict)
logging.debug("Saved raw data interp dictionary")

save(destination+'flux_interp_dict.hkl', flux_interp_dict)
logging.debug("Saved flux interp dictionary")
save(destination+'div_flux_dict.hkl', div_flux_dict)
logging.debug("Saved flux divergence interp")

save(destination+'flux_vert_interp_dict.hkl', flux_vert_interp_dict)
logging.debug("Saved flux vert dictionary")
save(destination+'div_vert_interp_dict.hkl', div_vert_interp_dict)
logging.debug("Saved flux divergence vert interp")


save(destination+"coord_dic.hkl", coord)

# import shutil
# for i in range(0,int(num)):
#  source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/reverse/'+dirc[1]+'_'+dirc[2]+str(i)+'/'
#  shutil.rmtree(source)

# logging.debug("Deleted individual 30 years !!")

logging.debug("Looks great !!")


