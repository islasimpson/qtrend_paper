import sys
import xarray as xr
import numpy as np


from qtrendutils import linfit_utils as linfit
from qtrendutils import calendar_utils as cal

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
topdir="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/qprofiles/"

dat = xr.open_dataset(topdir+'Q_sw_era5_plev.nc')
#dat = dat.groupby('time.year').mean('time')
dat = cal.calcannualmean(dat)
dat = dat.load()

dattrend = xr.apply_ufunc(linfit.compute_slope, dat, vectorize=True, 
            input_core_dims=[['year']])*dat.year.size

dattrend.to_netcdf(pathout+'qprofile_trend_ERA5.nc')
