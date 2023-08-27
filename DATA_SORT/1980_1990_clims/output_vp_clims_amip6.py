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

amip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/AMIP6/amip6csvinfo.csv')
models = amip6models['Model']
nmems = amip6models['Nmem']
nmemmax = np.max(amip6models['Nmem'])

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/AMIP6/vapor_pressures/"

amip6datt=[]

for index, modname in models.iteritems():
    print(modname)
    nmemsp = amip6models.loc[index, "Nmem"]
    dat = xr.open_dataset(datpath+'vaporpressures_'+modname+'.nc')
    #dat = dat.groupby('time.year').mean('time')
    dat = cal.calcannualmean(dat)
    dat = dat.sel(year=slice(1980,1990))
    dat = dat.mean('year')

    amip6datt.append(dat)

amip6dat = xr.concat(amip6datt, dim='model', coords='minimal', compat='override')
amip6dat = amip6dat.assign_coords(model=('model',models))

amip6dat.to_netcdf(pathout+'vpclims_AMIP6.nc')


