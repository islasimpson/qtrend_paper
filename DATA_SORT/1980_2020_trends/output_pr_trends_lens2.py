import sys
import xarray as xr
import numpy as np
import glob

from qtrendutils import lensread_utils as lensread
from qtrendutils import linfit_utils as linfit
from qtrendutils import calendar_utils as cal

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
topdir="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/precip/"

memnames = lensread.lens2memnamegen(100)
filelist = [ glob.glob(topdir+'*'+memstr+'*.nc')[0] for memstr in memnames ]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim='M')
#dat = dat.groupby('time.year').mean('time')
dat = cal.calcannualmean(dat)
dat = dat.load()

dattrend = xr.apply_ufunc(linfit.compute_slope, dat, vectorize=True,
               input_core_dims=[['year']])*dat.year.size

dattrend.to_netcdf(pathout+'prtrends_LENS2.nc')
