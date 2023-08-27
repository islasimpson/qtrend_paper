import importlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from math import nan
import scipy.stats as stats
import matplotlib as mpl

from CASutils import linfit_utils as linfit

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"

cmip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/CMIP6/cmip6csvinfo.csv')
models = cmip6models['Model']
nmems = cmip6models['Nmem']
nmemmax = np.max(cmip6models['Nmem'])

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP6/pr/"

cmip6datt=[]

for index, modname in models.iteritems():
    print(modname)
    nmemsp = cmip6models.loc[index, "Nmem"]
    dat = xr.open_dataset(datpath+'pr_'+modname+'.nc')
    dat = dat.groupby('time.year').mean('time')
    dat = dat.sel(year=slice(1980,2020))
    dattrend = xr.apply_ufunc(linfit.compute_slope, dat, vectorize=True,
                input_core_dims=[['year']])*dat.year.size

    cmip6datt.append(dattrend)

cmip6dat = xr.concat(cmip6datt, dim='model', coords='minimal', compat='override')
cmip6dat = cmip6dat.assign_coords(model=('model',models))

cmip6dat.to_netcdf(pathout+'prtrends_CMIP6.nc')


