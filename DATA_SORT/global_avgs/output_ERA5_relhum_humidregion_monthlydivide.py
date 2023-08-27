# First ordering into arid to humid months and then taking spatial averages over the
# most arid 6 months of the year and the most humid 6 months of the year
import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd
import sys

from qtrendutils import averaging_utils as avg
from qtrendutils import calendar_utils as cal

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/global_avgs/"

landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
oceanfrac = xr.where( np.isnan(landfrac), 1, nan)


pet_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020.nc")
pet_tc = pet_tc.sel(time=slice("1980-01-01","2020-12-31"))
pet_tc_monthly = pet_tc.groupby('time.month').mean('time')
pet_tc_am = cal.calcannualmean(pet_tc).mean('year')


pet_tc = pet_tc.mean('time')
ppt_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PPT_Terraclim_1980_2020.nc")
ppt_tc = ppt_tc.sel(time=slice("1980-01-01","2020-12-31"))
ppt_tc_monthly = ppt_tc.groupby('time.month').mean('time')
ppt_tc_am = cal.calcannualmean(ppt_tc).mean('year')

aridity_monthly = ppt_tc_monthly.ppt / pet_tc_monthly.pet
aridity_monthly['lon'] = landfrac.lon ; aridity_monthly['lat'] = landfrac.lat

aridity_am = ppt_tc_am / pet_tc_am
aridity_am['lon'] = landfrac.lon ; aridity_am['lat'] = landfrac.lat

aridity_am = aridity_am*landfrac
aridity_monthly = aridity_monthly*landfrac

aridmask = landfrac.where( (aridity_am > 0.05) & (aridity_am < 0.5), nan)
humidmask = landfrac.where( (aridity_am >= 0.5), nan)
hyperaridmask = landfrac.where( (aridity_am <= 0.05), nan)

aridmask = aridmask.where( aridmask.lat > -60, nan)
humidmask = humidmask.where( humidmask.lat > -60, nan)
hyperaridmask = hyperaridmask.where(hyperaridmask.lat > -60, nan)

landmask = landfrac.where( landfrac.lat > -60, nan)

# mask the monthly aridity indices to have limite the aridity index to 10 to prevent
# blow ups when dividing by small numbers
aridity_monthly = aridity_monthly.where(
      ( (aridity_monthly < 10) & (np.isfinite(aridity_monthly)) & ~np.isnan(aridity_monthly)), 10)


#---Figure out which data points to drop
# These will be ones where the climatological aridity is nan or infinite
# and where the trend in gpcc is nan
# also dropping all points that are not classed as humid
aridity_stack = aridity_am.stack(z=('lon','lat'))
gpcc_clim = xr.open_dataset(
      "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/"
       +"gpccclim_monthly.nc")
gpcc_clim = gpcc_clim.mean('month', skipna=False)
gpcc_clim['lon'] = landfrac.lon ; gpcc_clim['lat'] = landfrac.lat
gpcc_clim_stack = gpcc_clim.stack(z=('lon','lat'))

humidmask_stack = humidmask.stack(z=('lon','lat'))

idrop = np.argwhere( np.isnan(np.array(aridity_stack[:])) |
                     ~np.isfinite(np.array(aridity_stack[:])) |
                     np.isnan(np.array(gpcc_clim_stack.pr) ) | 
                     np.isnan(np.array(humidmask_stack)) )[:,0]



era5 = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/vaporpressures_ERA5.nc").sel(time=slice("1980-01","2020-12"))
era5 = era5.relhum
era5['lon'] = landfrac.lon ; era5['lat'] = landfrac.lat

# days in month weights for calculating averages
dayweights = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
dayweights = np.tile(dayweights, era5.lat.size*era5.lon.size)
dayweights = np.reshape(dayweights, [era5.lat.size, era5.lon.size,12])
dayweights = np.moveaxis(dayweights, 2, 0)
dayweights = xr.DataArray(dayweights,
              coords=[np.arange(0,12,1), era5.lat, era5.lon], 
              dims=['time','lat','lon'], name='dayweights')




era5_monthly = cal.group_monthly2yearly(era5)
era5_clim = era5_monthly.sel(year=slice(1980,1989)).mean('year')
era5_monthly =  (era5_monthly - era5_clim) 

aridity_monthly_stack = aridity_monthly.stack(loc=('lon','lat'))
era5_stack = era5_monthly.stack(loc=('lon','lat'))
dayweights_stack = dayweights.stack(loc=('lon','lat'))

aridity_monthly_stack = aridity_monthly_stack.drop_isel(loc=idrop)
era5_stack = era5_stack.drop_isel(loc=idrop)
dayweights_stack = dayweights_stack.drop_isel(loc=idrop)
lons = era5_stack.lon.rename('lons')
lats = era5_stack.lat.rename('lats')

aridity_monthly_stack = aridity_monthly_stack.reset_index('loc')
era5_stack = era5_stack.reset_index('loc')
dayweights_stack = dayweights_stack.reset_index('loc')
lons = lons.reset_index('loc')
lats = lats.reset_index('loc')

aridity_monthly_stack = aridity_monthly_stack.reset_coords( ('lon','lat'), drop=True)
era5_stack = era5_stack.reset_coords( ('lon','lat'), drop=True)
dayweights_stack = dayweights_stack.reset_coords( ('lon','lat'), drop=True)
lons = lons.reset_coords( ('lon','lat'), drop=True)
lats = lats.reset_coords( ('lon','lat'), drop=True)

dayweights_stack = np.take_along_axis(dayweights_stack.values, 
                       aridity_monthly_stack.argsort(axis=0), axis=0)

aridity_monthly_stack_expand = aridity_monthly_stack.expand_dims({'year':era5_stack.year.size},0)

era5_sorted = np.take_along_axis(era5_stack.values, aridity_monthly_stack_expand.argsort(axis=1), axis=1)


dayweights_stack_xr = xr.DataArray(dayweights_stack,
             coords=[era5_stack.time, np.arange(0,lons.size,1)],
             dims=['time','loc'], name='day_weights')


era5_sorted_xr = xr.DataArray(era5_sorted, 
                  coords=[era5_stack.year, era5_stack.time, np.arange(0,lons.size,1)],
                  dims=['year','time','loc'], name='era5_sorted')

dayweights_arid = dayweights_stack_xr.isel(time=slice(0,6))
dayweights_humid = dayweights_stack_xr.isel(time=slice(6,12))

era5_arid = era5_sorted_xr.isel(time=slice(0,6))
era5_humid = era5_sorted_xr.isel(time=slice(6,12))


#era5_aridm = era5_arid.mean('time')
#era5_humidm = era5_humid.mean('time')
#era5_am = era5_sorted_xr.mean('time')

era5_aridm = (era5_arid*dayweights_arid).sum('time') / dayweights_arid.sum('time')
era5_humidm = (era5_humid*dayweights_humid).sum('time') / dayweights_humid.sum('time')
era5_am = (era5_sorted_xr*dayweights_stack_xr).sum('time') / dayweights_stack_xr.sum('time')

def cosweight(darray, lats):
    weights = np.cos(np.deg2rad(lats))
    darrayw = darray.weighted(weights)
    darraym = darrayw.mean('loc', skipna=True)
    return darraym

era5_aridm_avg = cosweight(era5_aridm, lats) 
era5_humidm_avg = cosweight(era5_humidm, lats)
era5_am_avg = cosweight(era5_am, lats)

era5_aridm_avg = era5_aridm_avg.rename('arid')
era5_humidm_avg = era5_humidm_avg.rename('humid')
era5_am_avg = era5_am_avg.rename('am')

era5_aridm_avg.to_netcdf(pathout+'humidregions_arid_humid_halfyears_relhum.nc')
era5_humidm_avg.to_netcdf(pathout+'humidregions_arid_humid_halfyears_relhum.nc', mode='a')
era5_am_avg.to_netcdf(pathout+'humidregions_arid_humid_halfyears_relhum.nc', mode='a')



#dat.to_netcdf(pathout+'ERA5_q.nc')
