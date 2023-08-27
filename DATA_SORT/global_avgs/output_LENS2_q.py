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

basepath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/vaporpressures/"
lens2 = xr.open_mfdataset(basepath+"*.nc", combine='nested', concat_dim='M').q
lens2 = lens2.sel(time=slice("1980-01-01","2020-12-31"))
lens2['lon'] = landfrac.lon ; lens2['lat'] = landfrac.lat

lens2_all = avg.cosweightlonlat(lens2, 0, 360, -90, 90)
lens2_ocean = avg.cosweightlonlat(lens2*oceanfrac,0,360,-90,90)
lens2_land = avg.cosweightlonlat(lens2*landfrac, 0, 360, -90, 90)
lens2_arid = avg.cosweightlonlat(lens2*aridmask, 0, 360, -90, 90)
lens2_humid = avg.cosweightlonlat(lens2*humidmask, 0, 360, -90, 90)
lens2_hyperarid = avg.cosweightlonlat(lens2*hyperaridmask,0,360,-90,90)

lens2_all = lens2_all.rename('global')
lens2_ocean = lens2_ocean.rename('ocean')
lens2_land = lens2_land.rename('allland')
lens2_arid = lens2_arid.rename('arid')
lens2_humid = lens2_humid.rename('humid')
lens2_hyperarid = lens2_hyperarid.rename('hyperarid')

dat = xr.merge([lens2_all, lens2_ocean, lens2_land, lens2_arid, lens2_humid, lens2_hyperarid])

dat.to_netcdf(pathout+'LENS2_q.nc')
