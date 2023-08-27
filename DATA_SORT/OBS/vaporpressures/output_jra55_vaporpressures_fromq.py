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
qpath = "/project/mojave/observations/JRA55/day/QREFHT/"
pspath = "/project/mojave/observations/JRA55/day/PS/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"

tdat = read.read_sfc(tpath+"*.nc", "1979-01-01","2020-12-31")
qdat = read.read_sfc(qpath+"*.nc", "1979-01-01","2020-12-31")
psdat = read.read_sfc(pspath+"*.nc", "1979-01-01","2020-12-31")

monyearstr = xr.DataArray(tdat.indexes['time'].strftime('%Y-%m'), coords=tdat.time.coords, 
              name='monyearstr')
qdat_monthly = qdat.groupby(monyearstr).mean('time')
tdat_monthly = tdat.groupby(monyearstr).mean('time')
psdat_monthly = psdat.groupby(monyearstr).mean('time')
time = tdat.time.groupby(monyearstr).mean('time')
time = time.rename(monyearstr='time')
qdat_monthly = qdat_monthly.rename(monyearstr='time')
tdat_monthly = tdat_monthly.rename(monyearstr='time')
psdat_monthly = psdat_monthly.rename(monyearstr='time')
qdat_monthly['time'] = time
tdat_monthly['time'] = time
psdat_monthly['time'] = time


svp = qcalcs.calcsvp(tdat_monthly.T2m)
svp = svp.rename('svp')

vp = qcalcs.calcvpfromhuss(qdat_monthly.Q2m, psdat_monthly.PS)
vp = vp.rename('vp')

svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})

svp.to_netcdf(pathout+'vaporpressures_JRA55_fromq.nc')
vp.to_netcdf(pathout+'vaporpressures_JRA55_fromq.nc', mode='a')


