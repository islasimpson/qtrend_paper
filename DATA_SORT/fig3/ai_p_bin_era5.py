import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd
from CASutils import calendar_utils as cal

import sys

#------Read in the landfrac
landfrac = xr.open_dataset("../LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
landfrac = landfrac.where(landfrac.lat > -60, nan) # omitting Antarctica

#------Read in the Aridity Index from TerraClim regridded onto the CESM grid
ppt = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PPT_Terraclim_1980_2020.nc")
ppt = ppt.sel(time=slice("1980-01-01","2020-12-31"))
pet = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020.nc")
pet = pet.sel(time=slice("1980-01-01","2020-12-31"))

ppt_am = cal.calcannualmean(ppt).mean('year')
pet_am = cal.calcannualmean(pet).mean('year')

aridity = ppt_am / pet_am
aridity['lon'] = landfrac.lon ; aridity['lat'] = landfrac.lat
aridity = aridity*landfrac

trendpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"

#------Figure out which data points to drop
# These will be ones where the climatological aridity is nan or infinite
# and where the trend in gpcc is nan
aridity_stack = aridity.stack(z=('lon','lat'))
gpcc_trend = xr.open_dataset(trendpath+'prtrends_GPCC.nc')
gpcc_trend_stack = gpcc_trend.stack(z=('lon','lat'))

idrop = np.argwhere( np.isnan(np.array(aridity_stack[:])) |  
                     ~np.isfinite(np.array(aridity_stack[:])) | 
                     np.isnan(np.array(gpcc_trend_stack.pr) ) )[:,0]


#------1980 to 2020 trends
trendpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
era5dat = xr.open_dataset(trendpath+'vptrends_ERA5.nc')
era5dat['lon'] = landfrac.lon ; era5dat['lat'] = landfrac.lat
era5vp_map = era5dat.vp*landfrac
era5relhum_map = era5dat.relhum*landfrac
era5tas_map = era5dat.T2m*landfrac
era5q_map = era5dat.q*landfrac

gpcp = xr.open_dataset(trendpath+'prtrends_GPCP.nc')
gpcp['lon']  = landfrac.lon ; gpcp['lat'] = landfrac.lat
gpcp = gpcp.pr*landfrac

gpcc = xr.open_dataset(trendpath+'prtrends_GPCC.nc')
gpcc['lon']  = landfrac.lon ; gpcc['lat'] = landfrac.lat
gpcc = gpcc.pr*landfrac

cruts = xr.open_dataset(trendpath+'prtrends_CRUTS.nc')
cruts['lon']  = landfrac.lon ; cruts['lat'] = landfrac.lat
cruts = cruts.pr*landfrac

best = xr.open_dataset(trendpath+'tastrends_BEST.nc')
best['lon'] = landfrac.lon ; best['lat'] = landfrac.lat
best = best.tas*landfrac




#------1980 to 1990 climatologies
era5clim = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/"+
          "vpclims_ERA5.nc")
era5clim['lon'] = landfrac.lon ; era5clim['lat'] = landfrac.lat
era5clim_vp = era5clim.vp*landfrac
era5clim_q = era5clim.q*landfrac

#------Percent trends in vp and q
era5vp_map_pcnt = (era5vp_map / era5clim_vp)*100.
era5q_map_pcnt = (era5q_map / era5clim_q)*100.


#-----Stack in longitude and latitude
aridity_stack = aridity.stack(z=("lon","lat"))
era5vp_stack = era5vp_map.stack(z=("lon","lat"))
era5vp_pcnt_stack = era5vp_map_pcnt.stack(z=("lon","lat"))
era5q_pcnt_stack = era5q_map_pcnt.stack(z=("lon","lat"))
era5relhum_stack = era5relhum_map.stack(z=("lon","lat"))
era5tas_stack = era5tas_map.stack(z=("lon","lat"))
gpcp_stack = gpcp.stack(z=("lon","lat"))
gpcc_stack = gpcc.stack(z=("lon","lat"))
cruts_stack = cruts.stack(z=("lon","lat"))
best_stack = best.stack(z=("lon","lat"))

#idrop = np.argwhere( np.isnan(np.array(aridity_stack)) | ~np.isfinite(np.array(aridity_stack)))[:,0]
#idrop = np.argwhere( np.isnan(np.array(aridity_stack)) | ~np.isfinite(np.array(aridity_stack)) 
#                     | np.isnan(np.array(gpcc_stack)))[:,0]

aridity_stack = aridity_stack.drop_isel(z=idrop)
era5vp_stack = era5vp_stack.drop_isel(z=idrop)
era5vp_pcnt_stack = era5vp_pcnt_stack.drop_isel(z=idrop)
era5q_pcnt_stack = era5q_pcnt_stack.drop_isel(z=idrop)
era5relhum_stack = era5relhum_stack.drop_isel(z=idrop)
era5tas_stack = era5tas_stack.drop_isel(z=idrop)
gpcp_stack = gpcp_stack.drop_isel(z=idrop)
gpcc_stack = gpcc_stack.drop_isel(z=idrop)
cruts_stack = cruts_stack.drop_isel(z=idrop)
best_stack = best_stack.drop_isel(z=idrop)

nsplits_ai = 30 # 50 aridity bins
nsplits_pr = 15 # 25 precipitation bins

area_era5_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='area_gpcp')
aridity_ptiles = xr.DataArray(
             np.zeros([nsplits_ai])*nan,
             coords=[np.arange(0,nsplits_ai,1)],
             dims=['ai'], name='aridity_ptiles')
vp_era5_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='vpera5_gpcp')
vp_pcnt_era5_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='vpera5_pcnt_gpcp')
q_pcnt_era5_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='qera5_pcnt_gpcp')
tas_era5_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='tasera5_gpcp')
relhum_era5_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='relhumera5_gpcp')
pr_gpcp_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='prgpcp_gpcp')
best_gpcp = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='best_gpcp')


area_era5_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='area_gpcc')
#aridity_gpcc = xr.DataArray(
#             np.zeros([nsplits_ai, nsplits_pr])*nan,
#             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
#             dims=['ai','pr'], name='aridity_gpcc')
vp_era5_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='vpera5_gpcc')
vp_pcnt_era5_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='vpera5_pcnt_gpcc')
q_pcnt_era5_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='qera5_pcnt_gpcc')
tas_era5_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='tasera5_gpcc')
relhum_era5_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='relhumera5_gpcc')
pr_gpcc_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='prgpcc_gpcc')
best_gpcc = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='best_gpcc')


area_era5_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='area_cruts')
#aridity_cruts = xr.DataArray(
#             np.zeros([nsplits_ai, nsplits_pr])*nan,
#             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
#             dims=['ai','pr'], name='aridity_cruts')
vp_era5_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='vpera5_cruts')
vp_pcnt_era5_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='vpera5_pcnt_cruts')
q_pcnt_era5_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='qera5_pcnt_cruts')
tas_era5_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='tasera5_cruts')
relhum_era5_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='relhumera5_cruts')
pr_cruts_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='prcruts_cruts')
best_cruts = xr.DataArray(
             np.zeros([nsplits_ai, nsplits_pr])*nan,
             coords=[np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['ai','pr'], name='best_cruts')

#----Cosine weighted mean over data points
def meanw(dat):
    meanw = np.sum(dat*np.cos(np.deg2rad(dat.lat))) / np.sum(np.cos(np.deg2rad(dat.lat)))
    return meanw


#-----sort by aridity
indices = np.argsort(aridity_stack)
aridity_aisort = aridity_stack[np.array(indices)]
vp_aisort = era5vp_stack[np.array(indices)]
vp_pcnt_aisort = era5vp_pcnt_stack[np.array(indices)]
q_pcnt_aisort = era5q_pcnt_stack[np.array(indices)]
relhum_aisort = era5relhum_stack[np.array(indices)]
tas_aisort = era5tas_stack[np.array(indices)]
gpcp_aisort = gpcp_stack[np.array(indices)]
gpcc_aisort = gpcc_stack[np.array(indices)]
cruts_aisort = cruts_stack[np.array(indices)]
best_aisort = best_stack[np.array(indices)]

#-----define area weighting for each grid point
area = np.ones(len(indices))*np.cos(np.deg2rad(aridity_aisort.lat))

#-----split into 50 equal area bins according to aridity index
cum_area = area.cumsum()/area.sum()
idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_ai, endpoint=False)[1:])
area_chunks = np.split(area, idx)

vp_chunks = np.split(vp_aisort,idx)
aridity_chunks = np.split(aridity_aisort, idx)
vp_pcnt_chunks = np.split(vp_pcnt_aisort, idx)
q_pcnt_chunks = np.split(q_pcnt_aisort, idx)
relhum_chunks = np.split(relhum_aisort, idx)
tas_chunks = np.split(tas_aisort, idx)
gpcp_chunks = np.split(gpcp_aisort, idx)
gpcc_chunks = np.split(gpcc_aisort, idx)
cruts_chunks = np.split(cruts_aisort, idx)
best_chunks = np.split(best_aisort, idx)


#-----split into 25 equal area bins according to precipitation trend
for i in np.arange(0,nsplits_ai,1):
    area_bin = area_chunks[i]
    vp_bin = vp_chunks[i]
    aridity_bin = aridity_chunks[i]
    vp_pcnt_bin = vp_pcnt_chunks[i]
    q_pcnt_bin = q_pcnt_chunks[i]
    relhum_bin = relhum_chunks[i]
    tas_bin = tas_chunks[i]
    gpcp_bin = gpcp_chunks[i]
    gpcc_bin = gpcc_chunks[i]
    cruts_bin = cruts_chunks[i]
    best_bin = best_chunks[i]

#    aridity_ptiles[i] = aridity_bin.mean('z')
    aridity_ptiles[i] = meanw(aridity_bin)

    #-----sort based on GPCP
    indices = np.argsort(gpcp_bin)
    area_prsort = area_bin[np.array(indices)]
#    aridity_prsort = aridity_bin[np.array(indices)]
    vp_prsort = vp_bin[np.array(indices)]
    vp_pcnt_prsort = vp_pcnt_bin[np.array(indices)]
    q_pcnt_prsort = q_pcnt_bin[np.array(indices)]
    relhum_prsort = relhum_bin[np.array(indices)]
    tas_prsort = tas_bin[np.array(indices)]
    pr_prsort = gpcp_bin[np.array(indices)]
    best_prsort = best_bin[np.array(indices)]
   

    cum_area = area_prsort.cumsum()/area_prsort.sum()
    idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_pr, endpoint=False)[1:])
    area_prchunks = np.split(area_prsort, idx)
#    aridity_prchunks = np.split(aridity_prsort, idx)
    vp_prchunks = np.split(vp_prsort, idx)
    vp_pcnt_prchunks = np.split(vp_pcnt_prsort, idx)
    q_pcnt_prchunks = np.split(q_pcnt_prsort, idx)
    relhum_prchunks = np.split(relhum_prsort, idx)
    tas_prchunks = np.split(tas_prsort, idx)
    pr_prchunks = np.split(pr_prsort, idx)
    best_prchunks = np.split(best_prsort, idx)

    for j in np.arange(0,nsplits_pr,1):
        area_era5_gpcp[i,j] = area_prchunks[j].sum('z')
#        aridity_gpcp[i,j] = aridity_prchunks[j].sum('z')
#        vp_era5_gpcp[i,j] = vp_prchunks[j].mean('z')
#        vp_pcnt_era5_gpcp[i,j] = vp_pcnt_prchunks[j].mean('z')
#        q_pcnt_era5_gpcp[i,j] = q_pcnt_prchunks[j].mean('z')
#        relhum_era5_gpcp[i,j] = relhum_prchunks[j].mean('z')
#        tas_era5_gpcp[i,j] = tas_prchunks[j].mean('z')
#        pr_gpcp_gpcp[i,j] = pr_prchunks[j].mean('z')
#        best_gpcp[i,j] = best_prchunks[j].mean('z')
        vp_era5_gpcp[i,j] = meanw(vp_prchunks[j])
        vp_pcnt_era5_gpcp[i,j] = meanw(vp_pcnt_prchunks[j])
        q_pcnt_era5_gpcp[i,j] = meanw(q_pcnt_prchunks[j])
        relhum_era5_gpcp[i,j] = meanw(relhum_prchunks[j])
        tas_era5_gpcp[i,j] = meanw(tas_prchunks[j])
        pr_gpcp_gpcp[i,j] = meanw(pr_prchunks[j])
        best_gpcp[i,j] = meanw(best_prchunks[j])



    #---------------------------

    #-----sort based on GPCC
    indices = np.argsort(gpcc_bin)
    area_prsort = area_bin[np.array(indices)]
#    aridity_prsort = aridity_bin[np.array(indices)]
    vp_prsort = vp_bin[np.array(indices)]
    vp_pcnt_prsort = vp_pcnt_bin[np.array(indices)]
    q_pcnt_prsort = q_pcnt_bin[np.array(indices)]
    relhum_prsort = relhum_bin[np.array(indices)]
    tas_prsort = tas_bin[np.array(indices)]
    pr_prsort = gpcc_bin[np.array(indices)]
    best_prsort = best_bin[np.array(indices)]

    cum_area = area_prsort.cumsum()/area_prsort.sum()
    idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_pr, endpoint=False)[1:])
    area_prchunks = np.split(area_prsort, idx)
#    aridity_prchunks = np.split(aridity_prsort, idx)
    vp_prchunks = np.split(vp_prsort, idx)
    vp_pcnt_prchunks = np.split(vp_pcnt_prsort, idx)
    q_pcnt_prchunks = np.split(q_pcnt_prsort, idx)
    relhum_prchunks = np.split(relhum_prsort, idx)
    tas_prchunks = np.split(tas_prsort, idx)
    pr_prchunks = np.split(pr_prsort, idx)
    best_prchunks = np.split(best_prsort, idx)

    for j in np.arange(0,nsplits_pr,1):
        area_era5_gpcc[i,j] = area_prchunks[j].sum('z')
#        aridity_gpcc[i,j] = aridity_prchunks[j].mean('z')
#        vp_era5_gpcc[i,j] = vp_prchunks[j].mean('z')
#        vp_pcnt_era5_gpcc[i,j] = vp_pcnt_prchunks[j].mean('z')
#        q_pcnt_era5_gpcc[i,j] = q_pcnt_prchunks[j].mean('z')
#        relhum_era5_gpcc[i,j] = relhum_prchunks[j].mean('z')
#        tas_era5_gpcc[i,j] = tas_prchunks[j].mean('z')
#        pr_gpcc_gpcc[i,j] = pr_prchunks[j].mean('z')
#        best_gpcc[i,j] = best_prchunks[j].mean('z')

        vp_era5_gpcc[i,j] = meanw(vp_prchunks[j])
        vp_pcnt_era5_gpcc[i,j] = meanw(vp_pcnt_prchunks[j])
        q_pcnt_era5_gpcc[i,j] = meanw(q_pcnt_prchunks[j])
        relhum_era5_gpcc[i,j] = meanw(relhum_prchunks[j])
        tas_era5_gpcc[i,j] = meanw(tas_prchunks[j])
        pr_gpcc_gpcc[i,j] = meanw(pr_prchunks[j])
        best_gpcc[i,j] = meanw(best_prchunks[j])


    #---------------------------

    #-----sort based on CRUTS 
    indices = np.argsort(cruts_bin)
    area_prsort = area_bin[np.array(indices)]
#    aridity_prsort = aridity_bin[np.array(indices)]
    vp_prsort = vp_bin[np.array(indices)]
    vp_pcnt_prsort = vp_pcnt_bin[np.array(indices)]
    q_pcnt_prsort = q_pcnt_bin[np.array(indices)]
    relhum_prsort = relhum_bin[np.array(indices)]
    tas_prsort = tas_bin[np.array(indices)]
    pr_prsort = cruts_bin[np.array(indices)]
    best_prsort = best_bin[np.array(indices)]

    cum_area = area_prsort.cumsum()/area_prsort.sum()
    idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_pr, endpoint=False)[1:])
    area_prchunks = np.split(area_prsort, idx)
#    aridity_prchunks = np.split(aridity_prsort, idx)
    vp_prchunks = np.split(vp_prsort, idx)
    vp_pcnt_prchunks = np.split(vp_pcnt_prsort, idx)
    q_pcnt_prchunks = np.split(q_pcnt_prsort, idx)
    relhum_prchunks = np.split(relhum_prsort, idx)
    tas_prchunks = np.split(tas_prsort, idx)
    pr_prchunks = np.split(pr_prsort, idx)
    best_prchunks = np.split(best_prsort, idx)


    for j in np.arange(0,nsplits_pr,1):
        area_era5_cruts[i,j] = area_prchunks[j].sum('z')
#        aridity_cruts[i,j] = aridity_prchunks[j].mean('z')
#        vp_era5_cruts[i,j] = vp_prchunks[j].mean('z')
#        vp_pcnt_era5_cruts[i,j] = vp_pcnt_prchunks[j].mean('z')
#        q_pcnt_era5_cruts[i,j] = q_pcnt_prchunks[j].mean('z')
#        relhum_era5_cruts[i,j] = relhum_prchunks[j].mean('z')
#        tas_era5_cruts[i,j] = tas_prchunks[j].mean('z')
#        pr_cruts_cruts[i,j] = pr_prchunks[j].mean('z')
#        best_cruts[i,j] = best_prchunks[j].mean('z')

        vp_era5_cruts[i,j] = meanw(vp_prchunks[j])
        vp_pcnt_era5_cruts[i,j] = meanw(vp_pcnt_prchunks[j])
        q_pcnt_era5_cruts[i,j] = meanw(q_pcnt_prchunks[j])
        relhum_era5_cruts[i,j] = meanw(relhum_prchunks[j])
        tas_era5_cruts[i,j] = meanw(tas_prchunks[j])
        pr_cruts_cruts[i,j] = meanw(pr_prchunks[j])
        best_cruts[i,j] = meanw(best_prchunks[j])



    #---------------------------


dat = xr.merge([area_era5_gpcp,vp_era5_gpcp,vp_pcnt_era5_gpcp,q_pcnt_era5_gpcp, 
                  relhum_era5_gpcp,tas_era5_gpcp,pr_gpcp_gpcp,best_gpcp,
                area_era5_gpcc,vp_era5_gpcc,vp_pcnt_era5_gpcc,q_pcnt_era5_gpcc,
                  relhum_era5_gpcc,tas_era5_gpcc,pr_gpcc_gpcc,best_gpcc,
                area_era5_cruts,vp_era5_cruts,vp_pcnt_era5_cruts,q_pcnt_era5_cruts,
                  relhum_era5_cruts,tas_era5_cruts,pr_cruts_cruts,best_cruts,aridity_ptiles])
 
dat.to_netcdf("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/fig3/"+
   "obsbinned_ai"+str(nsplits_ai)+"_pr"+str(nsplits_pr)+"_weighted.nc")

















