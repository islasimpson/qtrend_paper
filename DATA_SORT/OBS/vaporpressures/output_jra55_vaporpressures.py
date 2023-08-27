#----Calculating JRA-55 vapor pressure from Dew Point Temperature
import xarray as xr
import numpy as np
import sys
import pandas as pd
import math
import xesmf as e

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs

tpath = "/project/mojave/observations/JRA55/day/TREFHT/"
tdpath = "/project/mojave/observations/JRA55/day/DEWP/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"

tdat = read.read_sfc(tpath+"*.nc", "1979-01-01","2020-12-31")
tddat_dep = read.read_sfc(tdpath+"*.nc", "1979-01-01","2020-12-31")
tddat = np.array(tdat.T2m) - np.array(tddat_dep.DEWP)
tddat = xr.DataArray(tddat, coords=tdat.T2m.coords, name='dewp')

monyearstr = xr.DataArray(tdat.indexes['time'].strftime('%Y-%m'), coords=tdat.time.coords, 
              name='monyearstr')
tddat_monthly = tddat.groupby(monyearstr).mean('time')
tdat_monthly = tdat.groupby(monyearstr).mean('time')
time = tdat.time.groupby(monyearstr).mean('time')
time = time.rename(monyearstr='time')
tddat_monthly = tddat_monthly.rename(monyearstr='time')
tdat_monthly = tdat_monthly.rename(monyearstr='time')
tddat_monthly['time'] = time
tdat_monthly['time'] = time

svp = qcalcs.calcsvp(tdat_monthly.T2m)
svp = svp.rename('svp')

vp = qcalcs.calcsvp(tddat_monthly)
vp = vp.rename('vp')

svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})

svp.to_netcdf(pathout+'vaporpressures_JRA55.nc')
vp.to_netcdf(pathout+'vaporpressures_JRA55.nc', mode='a')


