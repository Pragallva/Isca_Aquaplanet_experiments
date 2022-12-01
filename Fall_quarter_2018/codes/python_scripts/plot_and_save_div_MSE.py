import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import pylab as py
#import Extract_ncfile_save_fluxes_radiation
#import Reload_save_interpolated as svintp
import matplotlib.cm as cm
from reverse_cmap import rcmap

import sys
dirc=sys.argv

import os
import errno

BuRd = rcmap(cm.RdBu)

import hickle as hkl

import logging
log_directory='/project2/tas1/pragallva/Spring_quarter_2018/codes/shell_script/log/'+dirc[1]+'_'+dirc[2]
logging.basicConfig(filename=log_directory+'.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logging.debug('...........plot_and_save_div_MSE.py.............')

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
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

dirc=sys.argv

####################
#### smoothening ###
####################

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

source='/project2/tas1/pragallva/Spring_quarter_2018/post_process_data/'+dirc[1]+'_'+dirc[2]+'/'
            
make_sure_path_exists("/project2/tas1/pragallva/Spring_quarter_2018/figures/")
make_sure_path_exists("/project2/tas1/pragallva/Spring_quarter_2018/figures/MSE_fluxes/")

fig_dest="/project2/tas1/pragallva/Spring_quarter_2018/figures/MSE_fluxes/"

dirc=sys.argv
#title=dirc[1]+"  "+dirc[2].split("_")[0]+"m"+"  "+dirc[2].split("_")[1]
#title =dirc[1]+"  "+dirc[2]

second_part=dirc[2].split("_")
if len(second_part)>1 :
    title = dirc[1]+"  "+second_part[0]+"  "+second_part[1]
else :
    title = dirc[1]+"  "+second_part[0]

div=load(source+"div_flux_dict.hkl")
flux=load(source+"flux_interp_dict.hkl")

NE=div['TOA_d']+div['SFC_u']-div['dhdt']
NE1=div['SE']+div['MM']+div['TE']
latns=div['latn']
  
m=range(13)
py.rc('text', usetex=True)
py.rc('font', family='serif', serif='Palatino',weight='bold')

def add(y):
   y1 = np.append(y,y[:,0,np.newaxis],axis=1)
   return y1

def add2(y):
   y1 = np.append(y[-1],y)
   return y1


n=50
v = np.linspace(-180, 180, n, endpoint=True)

def plot_div() :
    
    fig=py.figure(figsize=(20, 10))
    py.suptitle(title+r" $Wm^{-2}$",fontsize=35,y=1.005)
           
    py.subplot(221)
    b = py.contourf(range(13),latns,add(NE), v, cmap=BuRd); 
    # py.colorbar()
    mticks = ['J','F','M','A','M','J','J','A','S','O','N','D']
    py.xticks((m), (mticks),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.ylim(-80,80)
    py.title(r'$F_{NE}$',fontsize=30)

    py.subplot(222)
    b = py.contourf(range(13),latns,add(div['MM']), v, cmap=BuRd); 
    # py.colorbar()
    mticks = ['J','F','M','A','M','J','J','A','S','O','N','D']
    py.xticks((m), (mticks),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.ylim(-80,80)
    py.title(r'$F_{MM}$',fontsize=30)

    py.subplot(223)
    b = py.contourf(range(13),latns,add(div['SE']), v, cmap=BuRd); 
    # py.colorbar()
    mticks = ['J','F','M','A','M','J','J','A','S','O','N','D']
    py.xticks((m), (mticks),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.ylim(-80,80)
    py.title(r'$F_{SE}$',fontsize=30)

    py.subplot(224)
    b = py.contourf(range(13),latns,add(div['TE']), v, cmap=BuRd); 
    # py.colorbar()
    py.ylim(-80,80)
    #py.plot(range(13),add2(zero_F_TE_N[:-1]),'ro-',lw=5.8,ms=10)
    mticks = ['J','F','M','A','M','J','J','A','S','O','N','D']
    py.xticks((m), (mticks),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.title(r'$F_{TE}$',fontsize=30)

    py.tight_layout()
    py.savefig(fig_dest+title+'_div.pdf')
    print fig_dest+title+'_div.pdf'
    #py.show()

plot_div()

NE=flux['TOA_d']+flux['SFC_u']-flux['dhdt']
NE1=flux['SE']+flux['MM']+flux['TE']
latns=flux['latn']
latnr=flux['latnr']

def plot_flux() :
    
    v = np.arange(-10, 11.0, 1.0)
    m=range(13)
    fig=py.figure(figsize=(20, 10))
    py.suptitle(title+" (PW)",fontsize=35,y=0.95)
           
    py.subplot(221)
    b = py.contourf(range(13),latnr,add(NE), v, cmap=BuRd); 
    # py.colorbar()
    c=py.contour(m,latnr,add(NE),v, colors='k',linewidths=1.0);
    py.clabel(c,  inline=1,fmt = '%1.1f',inline_spacing=40, fontsize=15)
    mticks = ['J','F','M','A','M','J','J','A','S','O','N','D','J']
    py.xticks((m), (mticks),fontsize=15)
    py.yticks(range(-80,81,20),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.ylim(-80,80)
    py.title('NE',fontsize=30)

    py.subplot(222)
    b = py.contourf(range(13),latns,add(flux['MM']), v, cmap=BuRd); 
    # py.colorbar()
    c=py.contour(m,latns,add(flux['MM']),v, colors='k',linewidths=1.0);
    py.clabel(c,  inline=1,fmt = '%1.1f',inline_spacing=40, fontsize=15)
    py.xticks((m), (mticks),fontsize=15)
    py.yticks(range(-80,81,20),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.ylim(-80,80)
    py.title('MM',fontsize=30)

    py.subplot(223)
    b = py.contourf(range(13),latns,add(flux['SE']), v, cmap=BuRd); 
    c=py.contour(m,latns,add(flux['SE']),v, colors='k',linewidths=1.0);
    py.clabel(c,  inline=1,fmt = '%1.1f',inline_spacing=40, fontsize=15)
    # py.colorbar()
    py.xticks((m), (mticks),fontsize=15)
    py.yticks(range(-80,81,20),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.ylim(-80,80)
    py.title('SE',fontsize=30)

    py.subplot(224)
    b = py.contourf(range(13),latns,add(flux['TE']), v, cmap=BuRd); 
    # py.colorbar()
    c=py.contour(m,latns,add(flux['TE']),v, colors='k',linewidths=1.0);
    py.clabel(c,  inline=1,fmt = '%1.1f',inline_spacing=40, fontsize=15)
    py.ylim(-80,80)
    #py.plot(range(13),add2(zero_F_TE_N[:-1]),'ro-',lw=5.8,ms=10)
    py.xticks((m), (mticks),fontsize=15)
    py.yticks(range(-80,81,20),fontsize=15)
    py.tick_params(labelsize=18,size=4,width=2)
    py.title('TE',fontsize=30)

    py.subplots_adjust(left=0.12, right=0.88, top=0.88, bottom=0.10, wspace=0.15, hspace=0.3)
    py.savefig(fig_dest+title+'_flux.pdf')
    print fig_dest+title+'_flux.pdf'
    #py.show()
    
plot_flux()
