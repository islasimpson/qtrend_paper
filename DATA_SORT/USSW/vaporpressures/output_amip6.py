import xarray as xr
import numpy as np
import pandas as pd
import sys
from math import nan

from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

amip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/AMIP6/amip6csvinfo.csv')
models = amip6models['Model']
nmems = amip6models['Nmem']
nmemmax = np.max(amip6models['Nmem'])

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/AMIP6/vapor_pressures/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/USSW/vaporpressures/"
timenew = pd.date_range("1980-01","2015-01", freq='M') 


amipdatt=[]
for index, modname in models.iteritems():
    print(modname)
    nmemsp = amip6models.loc[index, "Nmem"]
    dat = xr.open_dataset(datpath+"vaporpressures_"+modname+".nc")
    dat = dat.sel(time=slice("1980-01","2014-12"))
    dat['time'] = timenew
 
    shpfile="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/shp/"+\
            "gadm36_USA_1.shp"
    mask = shp.maskgen(shpfile, dat, 
           ['California','Nevada','Utah','Arizona','New Mexico','Colorado'])

    datmaskd = avg.cosweightlonlat(dat*mask,0,360,-90,90)
    amipdatt.append(datmaskd)

amipdat = xr.concat(amipdatt, dim='model', coords='minimal', compat='override')
amipdat = amipdat.assign_coords(model=("model",models))

amipdat.to_netcdf(pathout+'vaporpressures_AMIP6_USSW.nc')



