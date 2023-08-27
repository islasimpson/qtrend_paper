import importlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import nan


from qtrendutils import readdata_utils as read
from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

obs=['CRUTS','GPCC','GPCP']

for iobs in obs:

    pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/USSW/precip/"
    path="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/precip/"+iobs+".nc"
    dat = read.read_sfc(path,"1980-01-01","2020-12-31")
   
    shpfile="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/shp/gadm36_USA_1.shp" 
    #shpfile = "/project/cas/islas/shapefiles/usa/gadm36_USA_1.shp"
    mask = shp.maskgen(shpfile, dat, 
        ['California','Nevada','Utah','Arizona','New Mexico','Colorado'])
    
    dat = dat.pr
    datmaskd = xr.DataArray(dat*mask, coords = dat.coords)
    dat_sw = avg.cosweightlonlat(datmaskd, 0, 360, -90, 90)
    dat_sw = dat_sw.rename('pr')
    dat_sw.to_netcdf(pathout+iobs+'_USSW.nc')
