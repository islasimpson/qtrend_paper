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

dat = xr.open_mfdataset("/project/cas/islas/python_savs/qtrend_paper/ERA5_test_vp_temporally/*.nc")
dat = dat.sel(time=slice("1980-01-01","2020-12-31"))

#dat = dat.groupby('time.year').mean('time')
dat = cal.calcannualmean(dat)
dat = dat.sel(year=slice(1980,2020))
dat = dat.load()


dattrend = xr.apply_ufunc(linfit.compute_slope, dat, vectorize=True,
            input_core_dims=[['year']])*dat.year.size

dattrend.to_netcdf(pathout+'vptrends_ERA5_daily_hourly.nc')
