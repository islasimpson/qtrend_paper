import importlib
import pandas as pd
import xarray as xr
import numpy as np
import xesmf as xe
import sys
import warnings
from scipy import signal
import datetime

from qtrendutils import readdata_utils as read
from qtrendutils import calendar_utils as cal

warnings.filterwarnings('ignore')

# Interpolation grid
cesmdat = xr.open_dataset(
    '../../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/tas/"
wgtfile=pathout+'wgtfile.nc'

dat = xr.open_dataset("/project/mojave/observations/OBS-TAS/best.tas.185001-202112.nc")
tas = dat.tas
#tas = tas.sel(time=slice("1980-01-01","2020-12-31"))

regridder = xe.Regridder(tas, grid_out, 'bilinear', periodic=True, reuse_weights=False,
          filename=wgtfile)

tas_rg = regridder(tas)

tas_rg = tas_rg.rename('tas')
time = cal.YYYYMM2date(tas_rg.time)
tas_rg['time'] = time
tas_rg = tas_rg.sel(time=slice("1980-01-01","2020-12-31"))

tas_rg.to_netcdf(pathout+"BEST_1980_2020.nc")

