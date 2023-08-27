import xarray as xr
import numpy as np
from math import nan

from qtrendutils import mapplot_utils as mymaps
from qtrendutils import averaging_utils as avg
from qtrendutils import linfit_utils as linfit
from qtrendutils import calendar_utils as cal
from qtrendutils import colormap_utils as mycolors
import pandas as pd

import sys

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/bin_monthly/"

#----Land fraction
landfrac = xr.open_dataset("../LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
landfrac = landfrac.where( landfrac.lat > -60, nan)

#----Aridity zone masks
pet_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020.nc")
pet_tc = pet_tc.sel(time=slice("1980-01-01","2020-12-31"))
pet_tc_monthly = pet_tc.groupby('time.month').mean('time')
pet_tc_am = cal.calcannualmean(pet_tc).mean('year')

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

pet_tc_monthly['lon'] = landfrac.lon ; pet_tc_monthly['lat'] = landfrac.lat
petmask = xr.DataArray(np.zeros([pet_tc_monthly.lat.size, pet_tc_monthly.lon.size]),
                  coords=[pet_tc_monthly.lat, pet_tc_monthly.lon],
                  dims=['lat','lon'], name='petmask')
petmask[:,:]=1
petmask = petmask.where( pet_tc_monthly.min('month').pet != 0., nan)
petmask_stack = petmask.stack(z=('lon','lat'))


#----Figure out which data points to drop
# These will be ones where the climatological aridity is nan or infinite
# and where the trend in gpcc is nan
aridity_stack = aridity_am.stack(z=('lon','lat'))
gpcc_trend = xr.open_dataset(
       "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
       +"prtrends_GPCC.nc")
gpcc_trend_stack = gpcc_trend.stack(z=('lon','lat'))
idrop = np.argwhere( np.isnan(np.array(aridity_stack[:])) | 
                     ~np.isfinite(np.array(aridity_stack[:])) |  
                     np.isnan(np.array(gpcc_trend_stack.pr) ) | 
                     np.isnan(np.array(petmask_stack) ) )[:,0]

#----ERA5 specific humidity trend
datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
era5vp_map = xr.open_dataset(datpath+'vptrends_ERA5_monthly.nc')
era5vp_map['lon'] = landfrac.lon ; era5vp_map['lat'] = landfrac.lat
era5vp_map = era5vp_map*landfrac

era5tas_map = xr.open_dataset(datpath+'tastrends_ERA5_monthly.nc')
era5tas_map['lon'] = landfrac.lon ; era5tas_map['lat'] = landfrac.lat
era5tas_map = era5tas_map*landfrac

era5q_map = xr.open_dataset(datpath+'qtrends_ERA5_monthly.nc')
era5q_map['lon'] = landfrac.lon ; era5q_map['lat'] = landfrac.lat
era5q_map = era5q_map*landfrac

#---ERA5 climatologies (for % calculation)
era5clim = xr.open_dataset(
   "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/"+
   "vpclims_ERA5_monthly.nc")
era5clim['lon'] = landfrac.lon ; era5clim['lat'] = landfrac.lat
era5clim_vp = era5clim.vp*landfrac
era5clim_q = era5clim.q*landfrac
era5clim_tas = era5clim.T2m*landfrac


era5vp_map_pcnt = (era5vp_map / era5clim_vp)*100.
era5q_map_pcnt = (era5q_map / era5clim_q)*100.


#----GPCC precipitation trend
gpcc_map = xr.open_dataset(datpath+'prtrends_GPCC_monthly.nc')
gpcc_map['lon'] = landfrac.lon ; gpcc_map['lat'] = landfrac.lat
gpcc_map = gpcc_map*landfrac

gpcc_map_am = xr.open_dataset(
   "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"+
   "prtrends_GPCC.nc")



#----CMIP6 regression coefficients
regcoefs = xr.open_dataset(pathout+'cmip6_monthly_regcoefs.nc')


#----Predict vp based on the monthly regression in CMIP6
vppredict = regcoefs.aglobal + regcoefs.bglobal*gpcc_map
era5vptrenddif_gpcc = era5vp_map.vp - np.array(vppredict.pr)

#----Now doing the binning
aridity_monthly_stack = aridity_monthly.stack(loc=('lon','lat')).rename('aridity_monthly')
aridity_am_stack = aridity_am.stack(loc=('lon','lat')).rename('aridity_am')

gpcc_stack = gpcc_map.stack(loc=('lon','lat'))
gpcc_stack_am = gpcc_map_am.stack(loc=('lon','lat')) 

vppredict_stack = vppredict.stack(loc=('lon','lat')).rename({'pr':'vppredict'})
trenddif_stack = era5vptrenddif_gpcc.stack(loc=('lon','lat')).rename('vptrenddif')
vp_stack = era5vp_map.stack(loc=('lon','lat'))
pr_stack = gpcc_map.stack(loc=('lon','lat'))
vp_pcnt_stack = era5vp_map_pcnt.stack(loc=('lon','lat'))
q_pcnt_stack = era5q_map_pcnt.stack(loc=('lon','lat'))
tas_stack = era5tas_map.stack(loc=('lon','lat'))
vp_clim_stack = era5clim_vp.stack(loc=('lon','lat'))
q_clim_stack = era5clim_q.stack(loc=('lon','lat'))
tas_clim_stack = era5clim_tas.stack(loc=('lon','lat'))

aridity_am_stack = aridity_am_stack.reset_index('loc')
aridity_monthly_stack = aridity_monthly_stack.reset_index('loc')

aridity_monthly_stack = aridity_monthly_stack.drop_isel(loc=idrop)
aridity_am_stack = aridity_am_stack.drop_isel(loc=idrop)
vppredict_stack = vppredict_stack.drop_isel(loc=idrop)
trenddif_stack = trenddif_stack.drop_isel(loc=idrop)
vp_stack = vp_stack.drop_isel(loc=idrop)
pr_stack = pr_stack.drop_isel(loc=idrop)
vp_pcnt_stack = vp_pcnt_stack.drop_isel(loc=idrop)
q_pcnt_stack = q_pcnt_stack.drop_isel(loc=idrop)
tas_stack = tas_stack.drop_isel(loc=idrop)
vp_clim_stack = vp_clim_stack.drop_isel(loc=idrop)
q_clim_stack = q_clim_stack.drop_isel(loc=idrop)
tas_clim_stack = tas_clim_stack.drop_isel(loc=idrop)

# test for monthly aridity indices that are still nan's because pet=0, set them to the maximum finite ai
aridity_monthly_stack = aridity_monthly_stack.where( np.isfinite(aridity_monthly_stack), nan)
maxai = np.nanmax(aridity_monthly_stack)
aridity_monthly_stack = aridity_monthly_stack.where( 
           ~np.isnan(aridity_monthly_stack) & np.isfinite(aridity_monthly_stack), maxai)
                             


lons = vp_stack.lon.rename('lons')
lats = vp_stack.lat.rename('lats')

nsplits_ai = 30 # 30 aridity bins

area_out = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='area')
vppredict = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='vppredict')
trenddif = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='trenddif')
vp = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='vp')
vp_pcnt = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='vp_pcnt')
q_pcnt = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='q_pcnt')
tas = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='tas')

pr = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='pr')
aridity = xr.DataArray(np.zeros([nsplits_ai, 12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='aridity')
vp_clim = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='vp_clim')
q_clim = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='q_clim')
tas_clim = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='tas_clim')

aridity_am = xr.DataArray(np.zeros([nsplits_ai])*nan,
       coords=[np.arange(0,nsplits_ai,1)], dims=['ai'], name='aridity_am')



#----cosine weighted mean over data points
def meanw(dat, lats):
    meanw = np.sum(dat*np.cos(np.deg2rad(np.array(lats)))) / np.sum(np.cos(np.deg2rad(np.array(lats))))
    return meanw

#----Sort along the month axis according to aridity
aridity_monthly_sorted = np.take_along_axis(aridity_monthly_stack.values, aridity_monthly_stack.argsort(axis=0), axis=0)
aridity_monthly_sorted = xr.DataArray(aridity_monthly_sorted, coords=aridity_monthly_stack.coords, 
                 dims=aridity_monthly_stack.dims, name='aridity_monthly')

vppredict_sorted = np.take_along_axis(vppredict_stack.vppredict.values, 
                 aridity_monthly_stack.argsort(axis=0), axis=0)
vppredict_sorted = xr.DataArray(vppredict_sorted, coords=aridity_monthly_stack.coords, 
                 dims=aridity_monthly_stack.dims, name='vppredict')

trenddif_sorted = np.take_along_axis(trenddif_stack.values, 
                   aridity_monthly_stack.argsort(axis=0), axis=0)
trenddif_sorted = xr.DataArray(trenddif_sorted, coords=aridity_monthly_stack.coords, 
                   dims=aridity_monthly_stack.dims, name='trenddif')

vp_sorted = np.take_along_axis(vp_stack.vp.values, aridity_monthly_stack.argsort(axis=0), axis=0)
vp_sorted = xr.DataArray(vp_sorted, coords=aridity_monthly_stack.coords, 
                 dims=aridity_monthly_stack.dims, name='vp')

pr_sorted = np.take_along_axis(pr_stack.pr.values, aridity_monthly_stack.argsort(axis=0), axis=0)
pr_sorted = xr.DataArray(pr_sorted, coords=aridity_monthly_stack.coords, 
                dims=aridity_monthly_stack.dims, name='pr')

vp_pcnt_sorted = np.take_along_axis(vp_pcnt_stack.vp.values, aridity_monthly_stack.argsort(axis=0), axis=0)
vp_pcnt_sorted = xr.DataArray(vp_pcnt_sorted, coords=aridity_monthly_stack.coords,
                   dims=aridity_monthly_stack.dims, name='vp_pcnt')

q_pcnt_sorted = np.take_along_axis(q_pcnt_stack.q.values, aridity_monthly_stack.argsort(axis=0), axis=0)
q_pcnt_sorted = xr.DataArray(q_pcnt_sorted, coords=aridity_monthly_stack.coords,
                   dims=aridity_monthly_stack.dims, name='q_pcnt')

tas_sorted = np.take_along_axis(tas_stack.tas.values, aridity_monthly_stack.argsort(axis=0), axis=0)
tas_sorted = xr.DataArray(tas_sorted, coords=aridity_monthly_stack.coords,
                   dims=aridity_monthly_stack.dims, name='tas')

vp_clim_sorted = np.take_along_axis(vp_clim_stack.values, aridity_monthly_stack.argsort(axis=0), axis=0)
vp_clim_sorted = xr.DataArray(vp_clim_sorted, coords=aridity_monthly_stack.coords,
                   dims=aridity_monthly_stack.dims, name='vp_clim')

q_clim_sorted = np.take_along_axis(q_clim_stack.values, aridity_monthly_stack.argsort(axis=0), axis=0)
q_clim_sorted = xr.DataArray(q_clim_sorted, coords=aridity_monthly_stack.coords,
                   dims=aridity_monthly_stack.dims, name='q_clim')

tas_clim_sorted = np.take_along_axis(tas_clim_stack.values, aridity_monthly_stack.argsort(axis=0), axis=0)
tas_clim_sorted = xr.DataArray(tas_clim_sorted, coords=aridity_monthly_stack.coords,
                   dims=aridity_monthly_stack.dims, name='tas_clim')



lons = lons.reset_index('loc')
lats = lats.reset_index('loc')


alldat = xr.merge([aridity_monthly_sorted, vppredict_sorted, trenddif_sorted, vp_sorted, 
                   vp_pcnt_sorted, q_pcnt_sorted, tas_sorted, pr_sorted, aridity_am_stack, 
                   lons, lats, vp_clim_sorted, q_clim_sorted, tas_clim_sorted])
alldatsort = alldat.sortby(alldat.aridity_am)

area = np.ones(alldatsort.lats.size)*np.cos(np.deg2rad(alldatsort.lats))


# Split into equal area bins according to aridity index
cum_area = area.cumsum() / area.sum()
idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_ai, endpoint=False)[1:])

area_chunks = np.split(area,idx)
lons_chunks = np.split(alldatsort.lons, idx)
lats_chunks = np.split(alldatsort.lats, idx)
vp_chunks = np.split(alldatsort.vp, idx, axis=1)
vp_pcnt_chunks = np.split(alldatsort.vp_pcnt, idx, axis=1)
q_pcnt_chunks = np.split(alldatsort.q_pcnt, idx, axis=1)
vppredict_chunks = np.split(alldatsort.vppredict, idx, axis=1)
trenddif_chunks = np.split(alldatsort.trenddif, idx, axis=1)
pr_chunks = np.split(alldatsort.pr, idx, axis=1)
aridity_chunks = np.split(alldatsort.aridity_monthly, idx, axis=1)
tas_chunks = np.split(alldatsort.tas, idx, axis=1)
vp_clim_chunks = np.split(alldatsort.vp_clim, idx, axis=1)
q_clim_chunks = np.split(alldatsort.q_clim, idx, axis=1)
tas_clim_chunks = np.split(alldatsort.tas_clim, idx, axis=1)
aridity_am_chunks = np.split(alldatsort.aridity_am, idx, axis=0)



for i in np.arange(0,nsplits_ai,1):
    for j in np.arange(0,12,1):
        area_out[i,j] = area_chunks[i].sum('loc')
        vppredict[i,j] = meanw( np.array(vppredict_chunks[i])[j,:], np.array(lats_chunks[i]))
        trenddif[i,j] = meanw( np.array(trenddif_chunks[i])[j,:], np.array(lats_chunks[i]))
        vp[i,j] = meanw( np.array(vp_chunks[i])[j,:], np.array(lats_chunks[i]))
        vp_pcnt[i,j] = meanw( np.array(vp_pcnt_chunks[i])[j,:], np.array(lats_chunks[i]))
        q_pcnt[i,j] = meanw( np.array(q_pcnt_chunks[i])[j,:], np.array(lats_chunks[i]))
        tas[i,j] = meanw( np.array(tas_chunks[i])[j,:], np.array(lats_chunks[i]))
        pr[i,j] = meanw( np.array(pr_chunks[i])[j,:], np.array(lats_chunks[i]))
        aridity[i,j] = meanw( np.array(aridity_chunks[i])[j,:], np.array(lats_chunks[i]))
        vp_clim[i,j] = meanw( np.array(vp_clim_chunks[i])[j,:], np.array(lats_chunks[i]))
        q_clim[i,j] = meanw( np.array(q_clim_chunks[i])[j,:], np.array(lats_chunks[i]))
        tas_clim[i,j] = meanw( np.array(tas_clim_chunks[i])[j,:], np.array(lats_chunks[i]))

    aridity_am[i] = meanw( np.array(aridity_am_chunks[i])[:], np.array(lats_chunks[i]))


datout = xr.merge([area_out, vp, vp_pcnt, q_pcnt, tas, pr, aridity, vppredict, trenddif,
                  vp_clim, q_clim, aridity_am, tas_clim])
datout.to_netcdf(pathout+'obs_binned_monthly_onlypetnonzero.nc')






































