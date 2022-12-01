from netCDF4 import Dataset, num2date
import metpy.calc as mcalc
from metpy.units import units
import numpy as np
import netCDF4 as nc
import sys
import os.path
import pylab as py
import matplotlib.cm as cm
import numpy.ma as ma

import sys
import os
import errno
dirc=sys.argv

import logging


import time as ti
start = ti.time()

CO2      = str(dirc[1]) ## "HC0_la50m_oc50m"
num      = str(dirc[2]) ## represents an annual year of data
var1     = str(dirc[3])

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

def as_str(x):
   for k, v in list(locals().iteritems()):
       if v is x:
             x_as_str = k
   return x_as_str

log_directory = quarter+'/codes/shell_scripts/convert_isentrope/log/'
make_sure_path_exists( path   = log_directory )

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

if CO2=='None':
     CO2 = 'echr0001'
source  = '/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/'+CO2+'/'
exp_dir = 'ATM_dm_pl_'+CO2+'_10'
one_year = source+exp_dir+num+'.nc'

pres_dir = 'ATM_dm_ml_'+CO2+'_10' 
pres_one_year = source+pres_dir+num+'.nc'

logging.basicConfig( filename = log_directory+exp_dir+num+'_isen.log',level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

destination = quarter+'/exp_data/echam/data_in_isentropic_coord/'+exp_dir+'/'
make_sure_path_exists(destination);

data = nc.Dataset(one_year,'r')

logging.debug(str(num)+'...........interpolation_from_pres_to_isentrope_coord.py.............')

logging.debug(str(num)+'imported nc file')

def get_data(variable=var1, f = one_year):
	ncfile = f
	v_var  = nc.Dataset(ncfile,'r')
	data   = v_var.variables
	each_model_name = data[variable][:]
	v_var.close()
	return each_model_name

var_units = {'u':units.meter / (units.second),'v':units.meter / (units.second),\
	     'omega':units.pascals / (units.second),'geopoth':units.meter,\
	     't':units.K, 'q':units.dimensionless, 'deltheta_delt': units.K}


last     =   len(get_data('time'))

ps       = get_data('aps', pres_one_year)[0:last]
plevs    = get_data('lev')* units.pascals


lat       = get_data('lat')
lon       = get_data('lon')
time      = get_data('time')[0:last]

logging.debug(str(num)+'imported coordinates')

def surface_pressure_correction(y):
	mask = (plevs.magnitude[None,:,None,None]<ps[:, None, :,:])
	mask = mask.astype(float)
	mask[mask<1]= np.nan # 
	y_correct = y*mask
        return y_correct


temp     = surface_pressure_correction(get_data('t')[0:last,...])* units.K
logging.debug(str(num)+'imported temp')

kappa     = 0.286
theta     = (temp.magnitude)*(10**5/plevs[None,:,None,None])**(kappa)
dt        = 24*60*60
theta_dot = np.gradient(theta, dt, axis=0)* units.K/units.second

if (var1 != 'dtheta_dp'):
        VAR1     = surface_pressure_correction(get_data(var1)[0:last,...])* var_units[var1]
        logging.debug(str(num)+'imported '+var1)
else:
        VAR1 = np.gradient(theta, plevs, axis=1)* units.K/units.pascals
        logging.debug(str(num)+'calculated '+var1)

# Array of desired pressure levels
isen_levels = np.array([260,280,285,290,295,300,305,310,\
                        320,325,330,335,340,350,360,400,425,450,475,500]) * units.K

def interp_to_isentrope(field1=VAR1): 
        
        isen1_levels1 = mcalc.isentropic_interpolation(isen_levels, plevs, temp, field1, axis=1)
	logging.debug(str(num)+'interpolation done')
        return isen1_levels1[0].magnitude, isen1_levels1[1].magnitude ## last one is pressure

### Create a netcdf file ###

## create a dimension
netcdf_file = destination+exp_dir+num+'.nc'

if os.path.isfile(netcdf_file):
   	dataset   = nc.Dataset(netcdf_file,'a',format='NETCDF4_CLASSIC')

        logging.debug('Interpolating ------> '+var1+'<------ to isentropic coordinates')
        varp1 = interp_to_isentrope(VAR1)[-1]
        logging.debug('Calculation successful: shape '+str(np.shape(varp1)) )

else:   

	dataset   = nc.Dataset(netcdf_file,'w',format='NETCDF4_CLASSIC')

	logging.debug('Interpolating ------> '+var1+',temp and press <------ to isentropic coordinates')
	pres_p,  varp1 = interp_to_isentrope(VAR1)
	logging.debug('Calculation successful: shape '+str(np.shape(varp1)) )


	levels    = dataset.createDimension('level', len(isen_levels))
	lats      = dataset.createDimension('lat', len(lat))
	lons      = dataset.createDimension('lon', len(lon))
	times     = dataset.createDimension('time', None)
	logging.debug(str(num)+'Create dimension')   

	## Netcdf attribute
	dataset.title = exp_dir+"-"+num+" year"
	lats     = dataset.createVariable('lat', np.float32, ('lat',))
	lons     = dataset.createVariable('lon', np.float32, ('lon',))
	levels   = dataset.createVariable('level', np.float32, ('level',))
	times    = dataset.createVariable('time', np.float64, ('time',)) 
	logging.debug(str(num)+'Create coordinates')

	lats[:]    = lat
	lons[:]    = lon
	times[:]   = time
	levels[:]  = isen_levels.magnitude

        logging.debug('Interpolating ------> theta dot <------ to isentropic coordinates')
        deltheta_delt_p = interp_to_isentrope(theta_dot)[-1]
        logging.debug('Calculation successful: pres shape '+str(np.shape(deltheta_delt_p)) )

	PRESs          =  dataset.createVariable('pres',            np.float32, ('time','level','lat','lon'))
        Deltheta_Delts =  dataset.createVariable('deltheta_delt',   np.float32, ('time','level','lat','lon'))

	PRESs[:,:,:,:]          = pres_p
        logging.debug(str(num)+'Assigned P dataset' )
        Deltheta_Delts[:,:,:,:] = deltheta_delt_p
        logging.debug(str(num)+'Assigned theta_dot dataset' )

# Create the actual 4-d variable
    
VAR1s   = dataset.createVariable(var1,   np.float32, ('time','level','lat','lon'))
logging.debug(str(num)+'Create variables')

VAR1s[:,:,:,:]   = varp1

logging.debug(str(num)+'assigned '+var1+' dataset')

dataset.close()

end = ti.time()

logging.debug(str(num)+'Program ends')
logging.debug(str(num)+'Total time taken - '+str(end-start))
logging.debug(str(num)+'------------------------')
