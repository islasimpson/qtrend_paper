import xarray as xr
import numpy as np
import pandas as pd
import sys
from math import nan

from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/vaporpressures/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/USSW/vaporpressures/"

dat = xr.open_mfdataset(datpath+'*.nc', combine='nested', concat_dim='M')

dat = dat.sel(time=slice("1980-01","2020-12"))

shpfile="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/shp/"+\
            "gadm36_USA_1.shp"

mask = shp.maskgen(shpfile, dat,
           ['California','Nevada','Utah','Arizona','New Mexico','Colorado'])

datmaskd = avg.cosweightlonlat(dat*mask, 0, 360, -90, 90)

datmaskd.to_netcdf(pathout+'vaporpressures_LENS2_USSW.nc')
