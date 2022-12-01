import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pylab as py
#import Extract_ncfile_save_fluxes_radiation
#import Reload_save_interpolated as svintp
import matplotlib.cm as cm
import sys
import os
import errno
dirc=sys.argv
from reverse_cmap import rcmap

# BuRd = rcmap(cm.RdBu)

import hickle as hkl

import logging
log_directory='/project2/tas1/pragallva/Fall_quarter_2017/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.debug('...........plot_stream_wind_etc.py.............')

#source_dirc=svintp.source
#exp_dirc=["aqua_2m/","aqua_20m/","land_rec20m/"]

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
        #logging.debug('destination folder created !')
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            #logging.debug('destination folder exists already!')
            raise

dirc=sys.argv

####################
#### smoothening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

import netCDF4 as nc
filename="/project2/tas1/pragallva/Fall_quarter_2017/land_files/Ruth_full.nc"
data=nc.Dataset(filename,'r')
land_mask=data.variables['land_mask'][:]
zsurf=data.variables['zsurf'][:]
lat=data.variables['lat'][:]
lon=data.variables['lon'][:]


source='/project2/tas1/pragallva/Fall_quarter_2017/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'
            
make_sure_path_exists("/project2/tas1/pragallva/Fall_quarter_2017/figures/")
make_sure_path_exists("/project2/tas1/pragallva/Fall_quarter_2017/figures/stream_wind_ht/")

fig_dest="/project2/tas1/pragallva/Fall_quarter_2017/figures/stream_wind_ht/"

dirc=sys.argv
#title=dirc[1]+"  "+dirc[2].split("_")[0]+"m"+"  "+dirc[2].split("_")[1]
#title =dirc[1]+"  "+dirc[2]

second_part=dirc[2].split("_")
if len(second_part)>1 :
    title = dirc[1]+"_"+second_part[0]+"_"+second_part[1]
else :
    title = dirc[1]+"_"+second_part[0]
    
def aoy(x): # average over years
   # x=x.reshape(5,1440,x.shape[1],x.shape[2],x.shape[3])
    return x

u_v           =load(source+"u_v.hkl")
logging.debug('loaded uv')
ht_temp_sphum =load(source+"ht_temp_sphum.hkl")
logging.debug('loaded ht,temp,sphum')
stream_dic    =load(source+"stream_dic.hkl")
logging.debug('loaded stream')
coord        =load(source+"coord_dic.hkl")    ## save the interpolated plevs here using something like my_dict.update({'key 2': 'value 2', 'key 3': 'value 3'})
logging.debug('loaded coord')



lat = coord['lat']
lon = coord['lon']
#plev= coord_dic['interp_plev']
plev= np.array([0.5 ,10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0])
u   = aoy(u_v['u'])
v   = aoy(u_v['v'])
ht  = aoy(ht_temp_sphum['ht'])
stm = aoy(stream_dic['stream_func'])
        
m=range(13)
py.rc('text', usetex=True)
py.rc('font', family='serif', serif='Palatino',weight='bold')

import reverse_cmap as rc
my_cmap_r  = rc.rcmap(cm.RdBu)


## Calculate mass stream function ##
a=6371e3 ## Radius of earth in meters

import scipy.integrate as integrate
def integ( x) :    
    lo=np.deg2rad(coord['lon']) ;  
    la=np.deg2rad(coord['lat']) ;
    
    l1=59; l2=151;
    lons_arr= np.squeeze(np.where((lon>l1) & (lon<l2)))

    v_wind = x[:,:,:,np.squeeze(lons_arr)].mean(axis=-1)
    v_wind = v_wind*np.cos(la)[None,None,:] 
    x=v_wind    
    l= plev
    ax=1
    int_x = integrate.cumtrapz(x, l ,axis=ax, initial=None) #  (This is basically integration from TOA to surface)                                        
    return int_x
mass_stream=1*(2.0*np.pi*a/10.0)*integ(v) # for vertical integration

logging.debug('calculated mass stream function')

def plot_monsoon(x, h=0):
    mo=7
    m=30*4*(mo-1)+16*4 
    twenty_days=20*4
    H1=850-1;H2=850+1
    H_lev= np.where((plev>H1) & (plev<H2))    
    L1=150-1;L2=150+1
    L_lev= np.where((plev>L1) & (plev<L2))    
    post_monsoon_H= np.squeeze(x[m:m+twenty_days,H_lev,:,:])
    post_monsoon_L= np.squeeze(x[m:m+twenty_days,L_lev,:,:])
    
    co_H        = (post_monsoon_H).mean(axis=0)
    co_L        = (post_monsoon_L).mean(axis=0)

    if h :
        co_H    = co_H*10 # height in the units of m^2/s^-2
        co_L    = co_L*10 # height in the units of m^2/s^-2
        vv      = np.arange(np.min(co_H),np.max(co_H),200)
        vc      = np.arange(np.min(co_L),np.max(co_L),2000)
        print vc
        print "np.min(co_L), np.max(co_L)", np.min(co_L), np.max(co_L)
        print "np.min(co_H),np.max(co_H)", np.min(co_H),np.max(co_H)
        #vv     = np.arange(-3E7,+3.5E7,+0.5E7)
    else :
        vv      = np.arange(-3E7,+3.5E7,+0.5E7)
        vc      = np.arange(np.min(co_L),np.max(co_L),+2E7)

    b           = py.contourf(coord['lon'], coord['lat'], co_H, vv, cmap=my_cmap_r,  extend='both'); 
    cbar=py.colorbar(b)    
    cbar.ax.tick_params(labelsize=15)
    c           = py.contour (coord['lon'], coord['lat'], co_L, vc,  colors='k',linewidths=0.5 )     
    #py.clabel(c, fmt = '%1.0E', fontsize=10, inline=1)    
    d           = py.contour(coord['lon'], coord['lat'], land_mask, colors='lightgrey')
    py.tick_params(labelsize=18,size=4,width=2)
    py.xticks(range(0,361,60),fontsize=15)
    py.yticks(range(-90,91,30),fontsize=15)
    
    
def plot_vert_monsoon(x, post=1):       
    if post :
        mo=7
        m=120*(mo-1)+16*4
        twenty_days=20*4
    else :
        mo=4
        m=120*(mo-1)+00*4
        twenty_days=20*4
    l1=59; l2=151;
    lons_arr= np.squeeze(np.where((lon>l1) & (lon<l2)))
    
    wind        = x[:,:,:,np.squeeze(lons_arr)].mean(axis=-1)
    z_wind      = np.squeeze(wind[m:m+twenty_days,:,:]).mean(axis=0)
    ms          = mass_stream[m:m+twenty_days,:,:].mean(axis=0)
    vv          = np.linspace(-60,60,21,endpoint=True)
    b           = py.contourf(coord['lat'], plev, z_wind , vv, cmap=my_cmap_r,  extend='both'); 
    cbar=py.colorbar(b)    
    cbar.ax.tick_params(labelsize=15)
    c           = py.contour(coord['lat'], plev[1:], ms, 15, colors='k',linewidths=0.5 ) 
    #; py.clabel(c, fmt = '%1.0E', fontsize=10, inline=1)        
    py.tick_params(labelsize=18,size=4,width=2)
    py.xlim(-60,60)
    py.xticks(range(-60,61,20),fontsize=15)
    py.gca().invert_yaxis()
    

def plot_Ruth() :
    
    fig=py.figure(figsize=(20, 10))
    py.suptitle(" Land-real other diagnostics",fontsize=32,y=0.99)
           
    py.subplot(221)
    # plot zonal wind and MMC stream function pre-monsoon 
    plot_vert_monsoon(u, post=0)
    py.title('u wind & MMC (pre onset)',fontsize=28)
    logging.debug('plot pre onset')

    py.subplot(222)
    # plot zonal wind and MMC stream function post-monsoon 
    plot_vert_monsoon(u, post=1)
    py.title('u wind & MMC (post onset)',fontsize=28)
    logging.debug('plot post onset')

    py.subplot(223)
    plot_monsoon(ht, h=1)
    # plot geopotential post-monsoon onset
    py.title('geo at 850 hPa and 150 hPa',fontsize=28)
    logging.debug('plot geo')

    py.subplot(224)
    plot_monsoon(stm)
    # plot stream function post-monsoon onset
    py.title('stream at 850 hPa and 150 hPa',fontsize=28)
    logging.debug('plot stream')

    #py.tight_layout()
    py.savefig(fig_dest+title+'_Ruth.pdf')
    print fig_dest+title+'_Ruth.pdf'
    py.subplots_adjust(left=0.12, right=0.88, top=0.88, bottom=0.10, wspace=0.15, hspace=0.5)
    #py.show()
    
plot_Ruth()

logging.debug('COMPLETE !!')
