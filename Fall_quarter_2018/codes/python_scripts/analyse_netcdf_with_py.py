import numpy as np
import netCDF4 as nc
import sys
import os.path

dirc=sys.argv
path=os.getcwd()
#filename=os.path.join(path,dirc[1],"aqua_20m.nc")
#filename=dirc[0]+'/atmos_4xdaily.nc'
filename="/project2/tas1/pragallva/Fall_quarter_2017/codes/land.nc"
print filename, 'ok'
data=nc.Dataset(filename,'r')
#var=data.variables[dirc[1]]
#print nc.num2date(var[:],units=var.units, calendar=var.calendar)
#print data.variables[dirc[1]]
print data.variables.keys()
#print data.variables['height']
#print data.variables['phalf']
#print data.variables['pfull']
