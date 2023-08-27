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

dat = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"+
 "vaporpressures_ERA5.nc")

dat = cal.calcannualmean(dat)
#dat = dat.groupby('time.year').mean('time')
dat = dat.sel(year=slice(1980,1990))
dat = dat.mean('year')

dat.to_netcdf(pathout+'vpclims_ERA5.nc')
