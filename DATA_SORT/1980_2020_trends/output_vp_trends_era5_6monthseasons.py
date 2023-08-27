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

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"


#----aridity indices for picking out the most arid and least arid 6 months of the yera
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
landfrac = landfrac.where( landfrac.lat > -60, nan)


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

# mask the monthly aridity indices to have limite the aridity index to 10 to prevent
# blow ups when dividing by small numbers
aridity_monthly = aridity_monthly.where(
      ( (aridity_monthly < 10) & (np.isfinite(aridity_monthly)) & ~np.isnan(aridity_monthly)), 10)
#-------------------------------------------------------

dat = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"+
 "vaporpressures_ERA5.nc")
dat['lon'] = landfrac.lon ; dat['lat'] = landfrac.lat

dat_monthly = cal.group_monthly2yearly(dat)
dat_clim = dat_monthly.sel(year=slice(1980,1989)).mean('year')
dat_monthly = dat_monthly - dat_clim
dat_q_pcnt = (dat_monthly.q / dat_clim.q)*100.
dat_t = dat_monthly.T2m


dayweights = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
dayweights = np.tile(dayweights, dat.lat.size*dat.lon.size)
dayweights = np.reshape(dayweights, [dat.lat.size, dat.lon.size,12])
dayweights = np.moveaxis(dayweights, 2, 0)
dayweights = xr.DataArray(dayweights,
              coords=[np.arange(0,12,1), dat.lat, dat.lon],
              dims=['time','lat','lon'], name='dayweights')

aridity_monthly_expand = aridity_monthly.expand_dims({'year':dat_monthly.year.size},0)
dat_q_pcnt_sort = np.take_along_axis(dat_q_pcnt.values, aridity_monthly_expand.argsort(axis=1), axis=1)
dat_t_sort = np.take_along_axis(dat_t.values, aridity_monthly_expand.argsort(axis=1), axis=1)
dayweights_sort = np.take_along_axis(dayweights.values, aridity_monthly.argsort(axis=0), axis=0)

dat_q_pcnt_sort_xr = xr.DataArray(dat_q_pcnt_sort,
     coords=[dat_monthly.year, dat_monthly.time, dat_monthly.lat, dat_monthly.lon],
     dims=['year','month','lat','lon'], name='q_pcnt_sorted')
dat_t_sort_xr = xr.DataArray(dat_t_sort,
     coords=[dat_monthly.year, dat_monthly.time, dat_monthly.lat, dat_monthly.lon],
     dims=['year','month','lat','lon'], name='t_sorted')
dayweights_sort_xr = xr.DataArray(dayweights_sort,
     coords=[dat_monthly.time, dat_monthly.lat, dat_monthly.lon],
     dims=['month','lat','lon'], name='dayweights')

dat_q_arid = dat_q_pcnt_sort_xr.isel(month=slice(0,6))
dat_q_humid = dat_q_pcnt_sort_xr.isel(month=slice(6,12))

dat_t_arid = dat_t_sort_xr.isel(month=slice(0,6))
dat_t_humid = dat_t_sort_xr.isel(month=slice(6,12))

dayweights_arid = dayweights_sort_xr.isel(month=slice(0,6))
dayweights_humid = dayweights_sort_xr.isel(month=slice(6,12))

dat_q_aridm = (dat_q_arid*dayweights_arid).sum('month') / dayweights_arid.sum('month')
dat_q_humidm = (dat_q_humid*dayweights_humid).sum('month') / dayweights_humid.sum('month')

dat_t_aridm = (dat_t_arid*dayweights_arid).sum('month') / dayweights_arid.sum('month')
dat_t_humidm = (dat_t_humid*dayweights_humid).sum('month') / dayweights_humid.sum('month')


qtrend_arid = xr.apply_ufunc(linfit.compute_slope, dat_q_aridm, vectorize=True,
                  input_core_dims=[['year']])*dat_q_aridm.year.size
qtrend_humid = xr.apply_ufunc(linfit.compute_slope, dat_q_humidm, vectorize=True,
                  input_core_dims=[['year']])*dat_q_humidm.year.size

ttrend_arid = xr.apply_ufunc(linfit.compute_slope, dat_t_aridm, vectorize=True,
                  input_core_dims=[['year']])*dat_t_aridm.year.size
ttrend_humid = xr.apply_ufunc(linfit.compute_slope, dat_t_humidm, vectorize=True,
                  input_core_dims=[['year']])*dat_t_humidm.year.size


qtrend_arid = qtrend_arid.rename('q_arid')
qtrend_humid = qtrend_humid.rename('q_humid')
ttrend_arid = ttrend_arid.rename('t_arid')
ttrend_humid = ttrend_humid.rename('t_humid')

qtrend_arid.to_netcdf(pathout+'ERA5_vptrends_6monthseas_humid_arid.nc')
qtrend_humid.to_netcdf(pathout+'ERA5_vptrends_6monthseas_humid_arid.nc', mode='a')
ttrend_arid.to_netcdf(pathout+'ERA5_vptrends_6monthseas_humid_arid.nc', mode='a')
ttrend_humid.to_netcdf(pathout+'ERA5_vptrends_6monthseas_humid_arid.nc', mode='a')
