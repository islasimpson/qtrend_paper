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

dat = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"+
 "vaporpressures_ERA5.nc")
dat = dat.sel(time=slice("1980-01-01","2020-12-31")).T2m

year = dat.time.dt.year
nyears = dat.time.size/12
lat = dat.lat ; lon = dat.lon


dat_monthly = np.array(dat).reshape([int(nyears), 12, lat.size, lon.size])
year_monthly = np.array(year).reshape([int(nyears), 12])
yearvals = year_monthly[:,0]

dat_monthly = xr.DataArray(dat_monthly , dims=['year','month','lat','lon'], 
                  coords=[yearvals, np.arange(1,13,1), lat, lon], name='tas')

dat_monthly_trend = xr.apply_ufunc(linfit.compute_slope, dat_monthly, vectorize=True,
                      input_core_dims=[['year']])*dat_monthly.year.size

dat_monthly_trend.to_netcdf(pathout+'tastrends_ERA5_monthly.nc')
