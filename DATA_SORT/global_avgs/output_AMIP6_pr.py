import xarray as xr
import numpy as np
import matplotlib.pyplot as plt 
from math import nan
import pandas as pd

from qtrendutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/global_avgs/"

# Land fraction
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
oceanfrac = xr.where( np.isnan(landfrac), 1, nan)

# Aridity zone masks
pet_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020.nc")
pet_tc = pet_tc.mean('time')
ppt_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PPT_Terraclim_1980_2020.nc")
ppt_tc = ppt_tc.mean('time')
aridity = ppt_tc.ppt / pet_tc.pet
aridity['lon'] = landfrac.lon ; aridity['lat'] = landfrac.lat

aridmask = landfrac.where( (aridity > 0.05) & (aridity < 0.5), nan)
humidmask = landfrac.where( (aridity >= 0.5), nan)
hyperaridmask = landfrac.where( (aridity <= 0.05), nan)

aridmask = aridmask.where( aridmask.lat > -60, nan)
humidmask = humidmask.where( humidmask.lat > -60, nan)
hyperaridmask = hyperaridmask.where(hyperaridmask.lat > -60, nan)

landmask = landfrac.where( landfrac.lat > -60, nan)



cmip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/AMIP6/amip6csvinfo.csv')
models = cmip6models['Model']
nmems = cmip6models['Nmem']
nmemmax = np.max(cmip6models['Nmem'])

datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/AMIP6/PRECIP/"
timenew = pd.date_range("1980-01","2015-01", freq='M')

cmip6datt=[]
for index, modname in models.iteritems():
    print(modname)
    nmemsp = cmip6models.loc[index, "Nmem"]
    dat = xr.open_dataset(datpath+"precip_"+modname+".nc").pr
    dat = dat.sel(time=slice("1980-01","2014-12"))
    dat['time'] = timenew
    dat['lon'] = landmask.lon ; dat['lat'] = landmask.lat

    dat = dat.load()

    dat_ocean = avg.cosweightlonlat(dat*oceanfrac, 0, 360, -90, 90)
    dat_land = avg.cosweightlonlat(dat*landfrac, 0, 360, -90, 90)
    dat_arid = avg.cosweightlonlat(dat*aridmask, 0, 360, -90, 90)
    dat_humid = avg.cosweightlonlat(dat*humidmask, 0, 360, -90, 90)
    dat_hyperarid = avg.cosweightlonlat(dat*hyperaridmask,0,360,-90,90)

    dat_ocean = dat_ocean.rename('ocean')
    dat_land = dat_land.rename('allland')
    dat_arid = dat_arid.rename('arid')
    dat_humid = dat_humid.rename('humid')
    dat_hyperarid = dat_hyperarid.rename('hyperarid')

    dat = xr.merge([dat_ocean, dat_land, dat_arid, dat_humid, dat_hyperarid])

    cmip6datt.append(dat)

cmip6 = xr.concat(cmip6datt, dim='model', coords='minimal', compat='override')
cmip6 = cmip6.assign_coords(model=("model",models))

cmip6.to_netcdf(pathout+'AMIP6_pr.nc')



