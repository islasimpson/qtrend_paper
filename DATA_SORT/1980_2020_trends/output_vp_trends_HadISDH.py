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

dat = xr.open_dataset("/project/mojave/observations/HadISDH/gridded/"+
 "HadISDH.lande.4.4.0.2021f_FLATgridHOM5by5_anoms9120.nc")
dat = dat.e_abs

dat_seascyc = dat.groupby('time.month').mean('time', skipna=True)

dat_deseas = dat.groupby('time.month') - dat_seascyc

# add the time axis back into the seasonal cycle for computing annual averages
dat_seascyc = dat_seascyc.rename({'month':'time'})
dat_seascyc['time'] = dat.time.sel(time=slice("1980-01-01","1980-12-31"))

dat_seascyc_am = cal.calcannualmean(dat_seascyc)

# annual averaged deseasonalized anomalies
dat_deseas_am = cal.calcannualmean(dat_deseas, skipna=True)




#dat = dat.groupby('time.year').mean('time')
#dat = cal.calcannualmean(dat, skipna=True)
#dat = cal.calcannualmean(dat, skipna=False)
#dat = dat.sel(year=slice(1980,2020))
dat_deseas_am = dat_deseas_am.sel(year=slice(1980,2020))


dattrend = xr.apply_ufunc(linfit.compute_slope, dat_deseas_am, vectorize=True,
            input_core_dims=[['year']])*dat_deseas_am.year.size
dattrend = dattrend.rename('vp')

dattrend.to_netcdf(pathout+'vptrends_HadISDH.nc')
