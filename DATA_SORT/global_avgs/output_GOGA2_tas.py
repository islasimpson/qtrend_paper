import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd

from qtrendutils import averaging_utils as avg

import sys

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/global_avgs/"

landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
oceanfrac = xr.where( np.isnan(landfrac), 1, nan)

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

basepath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/GOGA2/vaporpressures/"
goga2 = xr.open_mfdataset(basepath+"*.nc", combine='nested', concat_dim='M').TREFHT
goga2 = goga2.sel(time=slice("1980-01-01","2020-12-31"))
goga2['lon'] = landfrac.lon ; goga2['lat'] = landfrac.lat

goga2_all = avg.cosweightlonlat(goga2, 0, 360, -90, 90)
goga2_ocean = avg.cosweightlonlat(goga2*oceanfrac,0,360,-90,90)
goga2_ocean_60sn = avg.cosweightlonlat(goga2*oceanfrac,0,360,-60,60)
goga2_land = avg.cosweightlonlat(goga2*landfrac, 0, 360, -90, 90)
goga2_arid = avg.cosweightlonlat(goga2*aridmask, 0, 360, -90, 90)
goga2_humid = avg.cosweightlonlat(goga2*humidmask, 0, 360, -90, 90)
goga2_hyperarid = avg.cosweightlonlat(goga2*hyperaridmask,0,360,-90,90)

goga2_all = goga2_all.rename('global')
goga2_ocean = goga2_ocean.rename('ocean')
goga2_ocean_60sn = goga2_ocean_60sn.rename('ocean_60sn')
goga2_land = goga2_land.rename('allland')
goga2_arid = goga2_arid.rename('arid')
goga2_humid = goga2_humid.rename('humid')
goga2_hyperarid = goga2_hyperarid.rename('hyperarid')

dat = xr.merge([goga2_all, goga2_ocean, goga2_ocean_60sn, goga2_land, goga2_arid, goga2_humid, goga2_hyperarid])

dat.to_netcdf(pathout+'GOGA2_tas.nc')
