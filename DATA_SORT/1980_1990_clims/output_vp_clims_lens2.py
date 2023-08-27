import importlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from math import nan
import scipy.stats as stats
import matplotlib as mpl
import glob

from qtrendutils import linfit_utils as linfit
from qtrendutils import calendar_utils as cal
from qtrendutils import lensread_utils as lensread

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/"
topdir="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/vaporpressures/"

memnames = lensread.lens2memnamegen(100)
filelist = [ glob.glob(topdir+'*'+memstr+'*.nc')[0] for memstr in memnames ]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim='M')
#dat = dat.groupby('time.year').mean('time')
dat = cal.calcannualmean(dat)

dat = dat.sel(year=slice(1980,1990))
dat = dat.mean('year')

dat.to_netcdf(pathout+'vpclims_LENS2.nc')

