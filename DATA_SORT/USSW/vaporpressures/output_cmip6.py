import xarray as xr
import numpy as np
import pandas as pd
import sys
from math import nan

from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

cmip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/CMIP6/cmip6csvinfo.csv')
models = cmip6models['Model']
nmems = cmip6models['Nmem']
nmemmax = np.max(cmip6models['Nmem'])

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP6/vapor_pressures/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/USSW/vaporpressures/"
timenew = pd.date_range("1980-01","2021-01", freq='M') 


cmipdatt=[]
for index, modname in models.iteritems():
    print(modname)
    nmemsp = cmip6models.loc[index, "Nmem"]
    dat = xr.open_dataset(datpath+"vaporpressures_"+modname+".nc")
    dat = dat.sel(time=slice("1980-01","2020-12"))
    dat['time'] = timenew
 
    shpfile="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/shp/"+\
            "gadm36_USA_1.shp"
    mask = shp.maskgen(shpfile, dat, 
           ['California','Nevada','Utah','Arizona','New Mexico','Colorado'])

    datmaskd = avg.cosweightlonlat(dat*mask,0,360,-90,90)
    cmipdatt.append(datmaskd)

cmipdat = xr.concat(cmipdatt, dim='model', coords='minimal', compat='override')
cmipdat = cmipdat.assign_coords(model=("model",models))

cmipdat.to_netcdf(pathout+'vaporpressures_CMIP6_USSW.nc')



