import importlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from math import nan
import scipy.stats as stats
import matplotlib as mpl

from qtrendutils import linfit_utils as linfit
from qtrendutils import calendar_utils as cal

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"

cmip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/CMIP6/cmip6csvinfo.csv')
models = cmip6models['Model']
nmems = cmip6models['Nmem']
nmemmax = np.max(cmip6models['Nmem'])

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP6/vapor_pressures/"

cmip6datt=[]

for index, modname in models.iteritems():
    print(modname)
    nmemsp = cmip6models.loc[index, "Nmem"]
    dat = xr.open_dataset(datpath+'vaporpressures_'+modname+'.nc')
    dat = dat.sel(time=slice("1980-01","2020-12")).tas

    year = dat.time.dt.year
    nyears = dat.time.size/12
    lat = dat.lat ; lon = dat.lon
    member = dat.member

    dat_monthly = np.array(dat).reshape([member.size,int(nyears),12,lat.size, lon.size])
    year_monthly = np.array(year).reshape([int(nyears),12])
    yearvals = year_monthly[:,0]

    dat_monthly = xr.DataArray(dat_monthly, dims=['member','year','month','lat','lon'],
                    coords=[member, yearvals, np.arange(1,13,1), lat, lon], name='tas')

    dat_monthly_trend = xr.apply_ufunc(linfit.compute_slope, dat_monthly, vectorize=True,
           input_core_dims=[['year']])*dat_monthly.year.size

    cmip6datt.append(dat_monthly_trend)

cmip6dat = xr.concat(cmip6datt, dim='model', coords='minimal', compat='override')
cmip6dat = cmip6dat.assign_coords(model=('model',models))

cmip6dat.to_netcdf(pathout+'tastrends_CMIP6_monthly.nc')


