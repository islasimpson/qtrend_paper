import importlib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from math import nan
import scipy.stats as stats
import matplotlib as mpl

from qtrendutils import linfit_utils as linfit
from qtrendutils import calendar_utils as cal


#----Get arid mask
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)

pet_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020.nc")
pet_tc = pet_tc.mean('time')
ppt_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PPT_Terraclim_1980_2020.nc")
ppt_tc = ppt_tc.mean('time')
aridity = ppt_tc.ppt / pet_tc.pet
aridity['lon'] = landfrac.lon ; aridity['lat'] = landfrac.lat

aridmask = landfrac.where( (aridity > 0.05) & (aridity < 0.5), nan)

aridmask = aridmask.where( aridmask.lat > -60, nan)
#-------------------


pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/rolling/"

dat = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"+
 "vaporpressures_ERA5.nc")
dat['lon'] = landfrac.lon ; dat['lat'] = landfrac.lat
dat = cal.calcannualmean(dat)
dat = dat.sel(year=slice(1980,2020))

dat = dat*aridmask
dat = dat.stack(z=['lon','lat']).dropna('z')
dat = dat.vp


yintrend=[19,21,23,25,27,29,31]
rollingtrend = xr.DataArray(np.zeros([len(yintrend), dat.year.size, dat.z.size]),
                 coords=[yintrend, dat.year, dat.z],
                 dims=['yintrend', 'year','z'], name='rolling_era5_trends')

for i in np.arange(0,len(yintrend),1):
    print(i)
    y = yintrend[i]
    rollingdat = dat.rolling(year=y, min_periods=y, center="True")
    rollingdat = rollingdat.construct("yinchunk")
    rollingtrendt = xr.apply_ufunc(linfit.compute_slope, rollingdat, vectorize=True,
      input_core_dims=[['yinchunk']])*rollingdat.yinchunk.size
    rollingtrend[i,:,:] = rollingtrendt 

#lons = dat.lon
#lats = dat.lat
rollingtrend = rollingtrend.reset_index('z')

rollingtrend.to_netcdf(pathout+'vptrends_ERA5_rolling.nc')
