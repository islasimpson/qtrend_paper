import importlib
import xarray as xr
import numpy as np
import xesmf as xe

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs
import sys


pathT = "/project/mojave/observations/ERA5/T2m/"
pathTD = "/project/mojave/observations/ERA5/TD2m/"
pathPS = "/project/mojave/observations/ERA5/PS/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"

datT = read.read_sfc(pathT+"*.nc", "1980-01","2020-12")
datTD = read.read_sfc(pathTD+"*.nc", "1980-01","2020-12")
datPS = read.read_sfc(pathPS+"*.nc", "1980-01","2020-12")

datPSclim = datPS.groupby('time.month').mean('time')
datPSclim_expand = datPSclim.expand_dims(dim={'year':np.arange(0,datPS.time.size/12.,1)})

datPS_stack = datPSclim_expand.stack(z=("year","month"))
datPS_stack = datPS_stack.transpose("z","lat","lon")
datPS_stack = datPS_stack.rename(z="time")
datPS_stack['time'] = datPS.time


#-----Regridding T and TD before the vp and svp calculation for consistency with CMIP
wgtfile=pathout+'wgtfile.nc'
cesmdat = xr.open_dataset(
     '../../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})
wgtreuse=False
regridder = xe.Regridder(datT.T2m, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                         filename=wgtfile)
datT_rg = regridder(datT)
datTD_rg = regridder(datTD)
datPS_rg = regridder(datPS_stack)

datT_rg = datT_rg.load()
datTD_rg = datTD_rg.load()
datPS_rg = datPS_rg.load()

print('after regridding')

#--Saturation vapor pressure
svp = qcalcs.calcsvp(datT_rg.T2m)
svp = svp.rename('svp')

#--Vapor pressure
vp = qcalcs.calcsvp(datTD_rg.TD2m)
vp = vp.rename('vp')

#--Vapor pressure deficit
vpd = svp - vp
vpd = vpd.rename('vpd')

#--Saturation specific humidity
sq = qcalcs.calcsq(datT_rg.T2m, datPS_rg.ps)

#--Specific humidity
q = qcalcs.calcsq(datTD_rg.TD2m, datPS_rg.ps)

#--Relative humidity
qrel = (q / sq)*100.
qrel = qrel.rename('relhum')

#--2m near surface air temperature
t2m = datT_rg.T2m
t2m = t2m.assign_attrs({'units':'K', 'long_name':'K'})

qrel = qrel.assign_attrs({'units':'%', 'long_name':'Relative Humidity'})
svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})
vpd = vpd.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure Deficit'})
q = q.assign_attrs({'units':'hPa', 'long_name':'Specific humidity'})
q = q.rename('q')

qrel.to_netcdf(pathout+'vaporpressures_ERA5_climps.nc')
svp.to_netcdf(pathout+'vaporpressures_ERA5_climps.nc', mode='a')
vp.to_netcdf(pathout+'vaporpressures_ERA5_climps.nc', mode='a')
vpd.to_netcdf(pathout+'vaporpressures_ERA5_climps.nc', mode='a')
t2m.to_netcdf(pathout+'vaporpressures_ERA5_climps.nc', mode='a')
q.to_netcdf(pathout+'vaporpressures_ERA5_climps.nc', mode='a') 
