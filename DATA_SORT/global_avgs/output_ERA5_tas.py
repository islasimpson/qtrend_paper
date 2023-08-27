import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd

from qtrendutils import averaging_utils as avg

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


era5 = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/vaporpressures_ERA5.nc").sel(time=slice("1980-01","2020-12")).T2m
era5['lon'] = landfrac.lon ; era5['lat'] = landfrac.lat

era5_all = avg.cosweightlonlat(era5, 0, 360, -90, 90)
era5_ocean = avg.cosweightlonlat(era5*oceanfrac,0,360,-90,90)
era5_ocean_60sn = avg.cosweightlonlat(era5*oceanfrac,0,360,-60,60)
era5_land = avg.cosweightlonlat(era5*landfrac,0,360,-90,90)
era5_arid = avg.cosweightlonlat(era5*aridmask,0,360,-90,90)
era5_humid = avg.cosweightlonlat(era5*humidmask,0,360,-90,90)
era5_hyperarid = avg.cosweightlonlat(era5*hyperaridmask,0,360,-90,90)

era5_all = era5_all.rename('global')
era5_ocean = era5_ocean.rename('ocean')
era5_ocean_60sn = era5_ocean_60sn.rename('ocean_60sn')
era5_land = era5_land.rename('allland')
era5_arid = era5_arid.rename('arid')
era5_humid = era5_humid.rename('humid')
era5_hyperarid = era5_hyperarid.rename('hyperarid')

dat = xr.merge([era5_all, era5_ocean, era5_ocean_60sn, era5_land, era5_arid, era5_humid, era5_hyperarid])

dat.to_netcdf(pathout+'ERA5_tas.nc')
