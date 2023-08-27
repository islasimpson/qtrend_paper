import importlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import nan

from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/USSW/vaporpressures/"
dat = xr.open_mfdataset("/project/cas/islas/python_savs/qtrend_paper/ERA5_test_vp_temporally/*.nc")
dat = dat.sel(time=slice("1980-01-01","2020-12-31"))
dat = dat.load()

shpfile="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/shp/gadm36_USA_1.shp"
mask = shp.maskgen(shpfile, dat,
   ['California','Nevada','Utah','Arizona','New Mexico','Colorado'])
datmaskd = dat*mask
dat_sw = avg.cosweightlonlat(datmaskd, 0, 360, -90, 90)
dat_sw.to_netcdf(pathout+'vaporpressures_ERA5_USSW_from_hourly_daily.nc')
