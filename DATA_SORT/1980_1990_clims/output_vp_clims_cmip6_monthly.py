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

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/"

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
    dat = dat.sel(time=slice("1980-01","1990-12"))
    dat = dat.groupby('time.month').mean('time')

    cmip6datt.append(dat)

cmip6dat = xr.concat(cmip6datt, dim='model', coords='minimal', compat='override')
cmip6dat = cmip6dat.assign_coords(model=('model',models))

cmip6dat.to_netcdf(pathout+'vpclims_CMIP6_monthly.nc')

