import importlib
import xarray as xr
import numpy as np
import xesmf as xe

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs

pathT = "/project/mojave/observations/ERA5/T2m/"
pathTD = "/project/mojave/observations/ERA5/TD2m/"
pathPS = "/project/mojave/observations/ERA5/PS/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"

datT = read.read_sfc(pathT+"*.nc", "1980-01","2021-12")
datTD = read.read_sfc(pathTD+"*.nc", "1980-01","2021-12")
datPS = read.read_sfc(pathPS+"*.nc", "1980-01","2021-12")

#-----Regridding T and TD before the vp and svp calculation for consistency with CMIP
#wgtfile=pathout+'wgtfile.nc'
#cesmdat = xr.open_dataset(
#     '../../cesm_grid.nc')
#grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})
#wgtreuse=False
#regridder = xe.Regridder(datT.T2m, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
#                         filename=wgtfile)
#datT_rg = regridder(datT)
#datTD_rg = regridder(datTD)
#datPS_rg = regridder(datPS)

#--Saturation vapor pressure
svp = qcalcs.calcsvp(datT.T2m)
svp = svp.rename('svp')

#--Vapor pressure
vp = qcalcs.calcsvp(datTD.TD2m)
vp = vp.rename('vp')

#--Vapor pressure deficit
vpd = svp - vp
vpd = vpd.rename('vpd')

#--Saturation specific humidity
sq = qcalcs.calcsq(datT.T2m, datPS.ps)

#--Specific humidity
q = qcalcs.calcsq(datTD.TD2m, datPS.ps)

#--Relative humidity
qrel = (q / sq)*100.
qrel = qrel.rename('relhum')

#--2m near surface air temperature
t2m = datT.T2m
t2m = t2m.assign_attrs({'units':'K', 'long_name':'K'})

qrel = qrel.assign_attrs({'units':'%', 'long_name':'Relative Humidity'})
svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})
vpd = vpd.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure Deficit'})

qrel.to_netcdf(pathout+'vaporpressures_ERA5_native.nc')
svp.to_netcdf(pathout+'vaporpressures_ERA5_native.nc', mode='a')
vp = vp.to_netcdf(pathout+'vaporpressures_ERA5_native.nc', mode='a')
vpd = vpd.to_netcdf(pathout+'vaporpressures_ERA5_native.nc', mode='a')
t2m = t2m.to_netcdf(pathout+'vaporpressures_ERA5_native.nc', mode='a')
