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

cmip6models=pd.read_csv('../CMIP6/cmip6csvinfo.csv')



#----Land fraction
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)

#----Aridity zone masks
pet_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020.nc")
pet_tc_monthly = pet_tc.groupby('time.month').mean('time')
pet_tc_am = cal.calcannualmean(pet_tc).mean('year')

ppt_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PPT_Terraclim_1980_2020.nc")
ppt_tc_monthly = ppt_tc.groupby('time.month').mean('time')
ppt_tc_am = cal.calcannualmean(ppt_tc).mean('year')

aridity_monthly = ppt_tc_monthly.ppt / pet_tc_monthly.pet
aridity_monthly['lon'] = landfrac.lon ; aridity_monthly['lat'] = landfrac.lat

aridity_am = ppt_tc_am / pet_tc_am
aridity_am['lon'] = landfrac.lon ; aridity_am['lat'] = landfrac.lat

landmask = landfrac.where( landfrac.lat > -60, nan)

#----ERA5 specific humidity trend
datpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
era5vp_map = xr.open_dataset(datpath+'vptrends_ERA5_monthly.nc')
era5vp_map['lon'] = landfrac.lon ; era5vp_map['lat'] = landfrac.lat

#----GPCC precipitation trend
gpcc_map = xr.open_dataset(datpath+'prtrends_GPCC_monthly.nc')
gpcc_map['lon'] = landfrac.lon ; gpcc_map['lat'] = landfrac.lat

#----CMIP6 data to work out regression coefficients as a function of month and location
cmip6vp_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/vptrends_CMIP6_monthly.nc')
cmip6vp_map['lon'] = landfrac.lon ; cmip6vp_map['lat'] = landfrac.lat

cmip6pr_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/prtrends_CMIP6_monthly.nc')*86400.
cmip6pr_map['lon'] = landfrac.lon ; cmip6pr_map['lat'] = landfrac.lat

cmip6vp_map_stack = cmip6vp_map.stack(z=('member','model'))
cmip6vp_map_stack = cmip6vp_map_stack.dropna('z')
cmip6vp_map_stack = cmip6vp_map_stack.rename({'vp':'vp_cmip6'})

cmip6pr_map_stack = cmip6pr_map.stack(z=('member','model'))
cmip6pr_map_stack = cmip6pr_map_stack.dropna('z')
cmip6pr_map_stack = cmip6pr_map_stack.rename({'pr':'pr_cmip6'})

a_global = np.zeros([cmip6pr_map.month.size, cmip6pr_map.lat.size, cmip6pr_map.lon.size])
b_global = np.zeros([cmip6pr_map.month.size, cmip6pr_map.lat.size, cmip6pr_map.lon.size])
a_global[:,:,:] = nan ; b_global[:,:,:] = nan

for imon in np.arange(0,12,1):
    print(imon)
    for ilon in np.arange(0,cmip6pr_map.lon.size,1):
        for ilat in np.arange(0,cmip6pr_map.lat.size,1):
            if (landfrac[ilat,ilon] == 1):
                a_t, b_t = linfit.linfit_xy(cmip6pr_map_stack.pr_cmip6.isel(month=imon, lon=ilon, lat=ilat), cmip6vp_map_stack.vp_cmip6.isel(month=imon, lon=ilon, lat=ilat))
                a_global[imon,ilat,ilon] = a_t ; b_global[imon,ilat,ilon] = b_t
                
a_global = xr.DataArray(a_global, coords=[cmip6pr_map_stack.month, cmip6pr_map_stack.lat, cmip6pr_map_stack.lon], dims=['month','lat','lon'], name='aglobal')
b_global = xr.DataArray(b_global, coords=[cmip6pr_map_stack.month, cmip6pr_map_stack.lat, cmip6pr_map_stack.lon], dims=['month','lat','lon'], name='bglobal')


#----Predict the vp based on the monthly regression in CMIP6 to derive the uncertainty range
vppredict = a_global +b_global*cmip6pr_map_stack.pr_cmip6
residuals = cmip6vp_map_stack.vp_cmip6 - vppredict
residualstdev = residuals.std(dim='z')
#vppredict = vppredict.rename({'vppredict':'vppredict_cmip6'})
vppredict = vppredict.rename('vppredict_cmip6')



#----Predict the ERA5 vp based on GPCC
vppredict_gpcc = a_global + b_global*gpcc_map
era5vptrenddif_gpcc = era5vp_map.vp - np.array(vppredict_gpcc.pr)

#----Now doing the binning
aridity_stack = aridity_monthly.stack(loc=('lon','lat')).rename('aridity_monthly')
aridity_am_stack = aridity_am.stack(loc=('lon','lat')).rename('aridity')

vppredict_cmip6_stack = vppredict.stack(loc=('lon','lat'))
vp_cmip6_stack = cmip6vp_map_stack.stack(loc=('lon','lat'))
pr_cmip6_stack = cmip6pr_map_stack.stack(loc=('lon','lat'))
vppredict_gpcc_stack = vppredict_gpcc.stack(loc=('lon','lat')).rename({'pr':'vppredict'})
trenddif_stack = era5vptrenddif_gpcc.stack(loc=('lon','lat')).rename('trenddif')
vp_stack = era5vp_map.stack(loc=('lon','lat'))
pr_stack = gpcc_map.stack(loc=('lon','lat'))

idrop = np.argwhere( np.isnan(np.array(aridity_stack[0,:])) | ~np.isfinite(np.array(aridity_stack[0,:])) | np.isnan(np.array(vppredict_gpcc_stack.vppredict[0,:])))[:,0]

vppredict_cmip6_stack = vppredict_cmip6_stack.drop_isel(loc=idrop)
vp_cmip6_stack = vp_cmip6_stack.drop_isel(loc=idrop)
pr_cmip6_stack = pr_cmip6_stack.drop_isel(loc=idrop)
aridity_stack = aridity_stack.drop_isel(loc=idrop)
aridity_am_stack = aridity_am_stack.drop_isel(loc=idrop)
trenddif_stack = trenddif_stack.drop_isel(loc=idrop)
vppredict_gpcc_stack = vppredict_gpcc_stack.drop_isel(loc=idrop)
vp_stack = vp_stack.drop_isel(loc=idrop)
pr_stack = pr_stack.drop_isel(loc=idrop)

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
pr = xr.DataArray(np.zeros([nsplits_ai,12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='pr')
aridity = xr.DataArray(np.zeros([nsplits_ai, 12])*nan,
       coords=[np.arange(0,nsplits_ai,1), np.arange(0,12,1)], dims=['ai','month'], name='aridity')

vp_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vpcmip6')

pr_cmip6 = xr.DataArray(
    np.zeros([cmip6pr_map.model.size, cmip6pr_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6pr_map.model, cmip6pr_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='prcmip6')

vppredict_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vppredict_cmip6')



#----cosine weighted mean over data points
def meanw(dat, lats):
    meanw = np.sum(dat*np.cos(np.deg2rad(np.array(lats)))) / np.sum(np.cos(np.deg2rad(np.array(lats))))
    return meanw


#----Sort along the month axis according to aridity

#--------------------OBS-----------------------------------

aridity_sorted = np.take_along_axis(aridity_stack.values, aridity_stack.argsort(axis=0), axis=0)
aridity_sorted = xr.DataArray(aridity_sorted, coords=aridity_stack.coords, dims=aridity_stack.dims, name='aridity_monthly')

vppredict_sorted = np.take_along_axis(vppredict_gpcc_stack.vppredict.values, aridity_stack.argsort(axis=0), axis=0)
vppredict_sorted = xr.DataArray(vppredict_sorted, coords=aridity_stack.coords, dims=aridity_stack.dims, name='vppredict')

trenddif_sorted = np.take_along_axis(trenddif_stack.values, aridity_stack.argsort(axis=0), axis=0)
trenddif_sorted = xr.DataArray(trenddif_sorted, coords=aridity_stack.coords, dims=aridity_stack.dims, name='trenddif')

vp_sorted = np.take_along_axis(vp_stack.vp.values, aridity_stack.argsort(axis=0), axis=0)
vp_sorted = xr.DataArray(vp_sorted, coords=aridity_stack.coords, dims=aridity_stack.dims, name='vp')

pr_sorted = np.take_along_axis(pr_stack.pr.values, aridity_stack.argsort(axis=0), axis=0)
pr_sorted = xr.DataArray(pr_sorted, coords=aridity_stack.coords, dims=aridity_stack.dims, name='pr')


alldat = xr.merge([aridity_sorted, vppredict_sorted, trenddif_sorted, vp_sorted, pr_sorted, aridity_am_stack, lons, lats])


alldatsort = alldat.sortby(alldat.aridity)

area = np.ones(alldat.lats.size)*np.cos(np.deg2rad(alldat.lats))

# Split into equal area bins according to aridity index
cum_area = area.cumsum() / area.sum()
idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_ai, endpoint=False)[1:])

area_chunks = np.split(area,idx)
lons_chunks = np.split(lons, idx)
lats_chunks = np.split(lats, idx)
vp_chunks = np.split(alldatsort.vp, idx, axis=1)
vppredict_chunks = np.split(alldatsort.vppredict, idx, axis=1)
trenddif_chunks = np.split(alldatsort.trenddif, idx, axis=1)
pr_chunks = np.split(alldatsort.pr, idx, axis=1)
aridity_chunks = np.split(alldatsort.aridity_monthly, idx, axis=1)

for i in np.arange(0,nsplits_ai,1):
    for j in np.arange(0,12,1):
        area_out[i,j] = area_chunks[i].sum('loc')
        vppredict[i,j] = meanw( np.array(vppredict_chunks[i])[j,:], np.array(lats_chunks[i]))
        trenddif[i,j] = meanw( np.array(trenddif_chunks[i])[j,:], np.array(lats_chunks[i]))
        vp[i,j] = meanw( np.array(vp_chunks[i])[j,:], np.array(lats_chunks[i]))
        pr[i,j] = meanw( np.array(pr_chunks[i])[j,:], np.array(lats_chunks[i]))
        aridity[i,j] = meanw( np.array(aridity_chunks[i])[j,:], np.array(lats_chunks[i]))


#MODEL
#
#vp_cmip6 = xr.DataArray(
#    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
#    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
#    dims=['model','member','ai','month'], name='vpcmip6')
#
#vppredict_cmip6 = xr.DataArray(
#    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
#    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
#    dims=['model','member','ai','month'], name='vppredictcmip6')
#
#pr_cmip6 = xr.DataArray(
#    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
#    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
#    dims=['model','member','ai','month'], name='prcmip6')
#
#aridity_check = xr.DataArray(
#    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
#    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0, nsplits_ai, 1), np.arange(0,12,1)],
#    dims=['model','member','ai','month'], name='aridity_check')
#
#
#
#
#for imodel in np.arange(0,cmip6vp_map.model.size,1):
#    print(imodel)
#    nmems = cmip6models.loc[imodel,"Nmem"]
#
#    for imem in np.arange(0,nmems,1):
#        vpdat = vp_cmip6_stack.isel(model=imodel, member=imem)
#        vppredictdat = vppredict_cmip6_stack.isel(model=imodel, member=imem)
#        prdat = pr_cmip6_stack.isel(model=imodel, member=imem)
#
#        aridity_sorted = np.take_along_axis(aridity_stack.values, aridity_stack.argsort(axis=0), axis=0)
#        aridity_sorted = xr.DataArray(aridity_sorted, coords=aridity_stack.coords, dims=aridity_stack.dims, name='aridity_monthly')
#
#        vp_sorted = np.take_along_axis(vpdat.values, aridity_stack.argsort(axis=0), axis=0)
#        vp_sorted = xr.DataArray(vp_sorted, coords=vpdat.coords, dims=vpdat.dims, name='vp_cmip6')
#
#        vppredict_sorted = np.take_along_axis(vppredictdat.values, aridity_stack.argsort(axis=0), axis=0)
#        vppredict_sorted = xr.DataArray(vppredict_sorted, coords=vppredictdat.coords, dims=vppredictdat.dims, name='vppredict_cmip6')
#
#        pr_sorted = np.take_along_axis(prdat.values, aridity_stack.argsort(axis=0), axis=0)
#        pr_sorted = xr.DataArray(pr_sorted, coords=prdat.coords, dims=prdat.dims, name='pr_cmip6')
#
#        alldat = xr.merge([vp_sorted, vppredict_sorted, pr_sorted, aridity_sorted, lons, lats, aridity_am_stack])
#
#        #----Sort by the annual mean aridity
#        alldatsort = alldat.sortby(alldat.aridity)
#
#        area = np.ones(alldatsort.lat.size)*np.cos(np.deg2rad(alldatsort.lat))
#
#        #----Split into equal area bins according to aridity index
#        cum_area = area.cumsum()/area.sum()
#        idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_ai, endpoint=False)[1:])
#
#        lons_chunks = np.split(lons, idx)
#        lats_chunks = np.split(lats, idx)
#        vp_chunks = np.split(alldatsort.vp_cmip6, idx, axis=1)
#        vppredict_chunks = np.split(alldatsort.vppredict_cmip6, idx, axis=1)
#        pr_chunks = np.split(alldatsort.pr_cmip6, idx, axis=1)
#        aridity_chunks = np.split(alldatsort.aridity_monthly, idx, axis=1)
#
#        for i in np.arange(0,nsplits_ai,1):
#            for j in np.arange(0,12,1):
#                vp_cmip6[imodel, imem, i, j] = meanw( np.array(vp_chunks[i])[j,:], np.array(lats_chunks[i]) ) 
#                vppredict_cmip6[imodel, imem, i, j] = meanw( np.array(vppredict_chunks[i])[j,:], np.array(lats_chunks[i]) ) 
#                pr_cmip6[imodel, imem, i, j] = meanw( np.array(pr_chunks[i])[j,:], np.array(lats_chunks[i]) )
#                aridity_check[imodel, imem, i, j] = meanw( np.array(aridity_chunks[i])[j,:], np.array(lats_chunks[i]) )
#

datout = xr.merge([area_out, vp, pr, aridity, vppredict, trenddif])
datout.to_netcdf(pathout+'obs_binned_monthly_OLD.nc')





