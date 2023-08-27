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

warnings.filterwarnings('ignore')

# Interpolation grid
cesmdat = xr.open_dataset(
    '../../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/precip/"
wgtfile=pathout+'wgtfile.nc'

# ---GPCP  (doesn't go back past 1979)
ybeg=1979 ; monbeg=1 ; yend=2021 ; monend=12 # dates for past period
nmonths = (yend-ybeg-1)*12 + (12-monbeg+1) + monend
fpath="/project/mojave/observations/OBS-PR/GPCP/gpcp.mon.mean.197901-202112.nc"
datebeg = str(ybeg)+"-"+str(monbeg).zfill(2)
dateend = str(yend)+"-"+str(monend).zfill(2)

pr = read.read_sfc(fpath, datebeg, dateend)

regridder = xe.Regridder(pr, grid_out, 'bilinear', periodic=True, reuse_weights=False,
                         filename=wgtfile)
pr_rg = regridder(pr)
pr_rg = pr_rg.rename({'precip':'pr'})
pr_rg.to_netcdf(pathout+"GPCP.nc")
