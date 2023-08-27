import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
import pandas as pd
from qtrendutils import calendar_utils as cal

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

#ppt = ppt.mean('time')
#pet = pet.mean('time')
ppt_am = cal.calcannualmean(ppt).mean('year')
pet_am = cal.calcannualmean(pet).mean('year')
aridity = ppt_am / pet_am
#aridity = ppt.ppt/pet.pet
aridity['lon'] = landfrac.lon ; aridity['lat'] = landfrac.lat
aridity = aridity*landfrac

#-------read in the CMIP6 models
cmip6models=pd.read_csv('../CMIP6/cmip6csvinfo.csv')



trendpath="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/"
cmip6dat = xr.open_dataset(trendpath+'vptrends_CMIP6.nc')
cmip6dat['lon'] = landfrac.lon ; cmip6dat['lat'] = landfrac.lat
cmip6vp_map = cmip6dat.vp*landfrac
cmip6relhum_map = cmip6dat.relhum*landfrac
cmip6tas_map = cmip6dat.tas*landfrac
cmip6q_map = cmip6dat.q*landfrac

cmip6pr_map = xr.open_dataset(trendpath+'prtrends_CMIP6.nc')*86400.
cmip6pr_map['lon'] = landfrac.lon ; cmip6pr_map['lat'] = landfrac.lat
cmip6pr_map = cmip6pr_map*landfrac
cmip6pr_map = cmip6pr_map.pr

clim_1980_1990 = xr.open_dataset(
             "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/vpclims_CMIP6.nc")
clim_1980_1990['lon'] = landfrac.lon ; clim_1980_1990['lat'] = landfrac.lat
cmip6vpclim_1980_1990 = clim_1980_1990.vp*landfrac
cmip6qclim_1980_1990 = clim_1980_1990.q*landfrac

cmip6vp_map_pcnt = (cmip6vp_map/cmip6vpclim_1980_1990)*100.
cmip6q_map_pcnt = (cmip6q_map/cmip6qclim_1980_1990)*100.

#-------Stack in longitude and latitude
aridity_stack = aridity.stack(z=("lon","lat"))
cmip6vp_stack = cmip6vp_map.stack(z=("lon","lat"))
cmip6relhum_stack = cmip6relhum_map.stack(z=("lon","lat"))
cmip6tas_stack = cmip6tas_map.stack(z=("lon","lat"))
cmip6pr_stack = cmip6pr_map.stack(z=("lon","lat"))
cmip6vp_pcnt_stack = cmip6vp_map_pcnt.stack(z=("lon","lat"))
cmip6q_pcnt_stack = cmip6q_map_pcnt.stack(z=("lon","lat"))

idrop = np.argwhere( np.isnan(np.array(aridity_stack)) | ~np.isfinite(np.array(aridity_stack)))[:,0]

aridity_stack = aridity_stack.drop_isel(z=idrop)
cmip6vp_stack = cmip6vp_stack.drop_isel(z=idrop)
cmip6relhum_stack = cmip6relhum_stack.drop_isel(z=idrop)
cmip6tas_stack = cmip6tas_stack.drop_isel(z=idrop)
cmip6pr_stack = cmip6pr_stack.drop_isel(z=idrop)
cmip6vp_pcnt_stack = cmip6vp_pcnt_stack.drop_isel(z=idrop)
cmip6q_pcnt_stack = cmip6q_pcnt_stack.drop_isel(z=idrop)

nsplits_ai = 30 # 30 aridity bins
nsplits_pr = 15 # 15 precipitation bins

area_cmip6 = xr.DataArray(
             np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size,nsplits_ai, nsplits_pr])*nan, 
             coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='area')
aridity_cmip6 = xr.DataArray(
             np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size,nsplits_ai, nsplits_pr])*nan,
             coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='aridity')
vp_cmip6 = xr.DataArray(
             np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size,nsplits_ai, nsplits_pr])*nan,
             coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='vpcmip6')
vp_cmip6_pcnt = xr.DataArray(
             np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size,nsplits_ai, nsplits_pr])*nan,
             coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='vpcmip6_pcnt')
q_cmip6_pcnt = xr.DataArray(
             np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size,nsplits_ai, nsplits_pr])*nan,
             coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='qcmip6_pcnt')
relhum_cmip6 = xr.DataArray(
             np.zeros([cmip6relhum_map.model.size, cmip6relhum_map.member.size,nsplits_ai, nsplits_pr])*nan,
             coords=[cmip6relhum_map.model, cmip6relhum_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='relhumcmip6')

tas_cmip6 = xr.DataArray(
          np.zeros([cmip6tas_map.model.size, cmip6tas_map.member.size, nsplits_ai, nsplits_pr])*nan,
          coords=[cmip6tas_map.model, cmip6tas_map.member, np.arange(0,nsplits_ai,1), np.arange(0,nsplits_pr,1)],
          dims=['model','member','ai','pr'], name='tascmip6')

pr_cmip6 = xr.DataArray(
             np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size,nsplits_ai, nsplits_pr])*nan,
             coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1),np.arange(0,nsplits_pr,1)],
             dims=['model','member','ai','pr'], name='prcmip6')

#----Cosine weighted mean over data points
def meanw(dat):
    meanw = np.sum(dat*np.cos(np.deg2rad(dat.lat))) / np.sum(np.cos(np.deg2rad(dat.lat)))
    return meanw


#-------Loop over models and members.
for imodel in np.arange(0,cmip6vp_stack.model.size,1):
    print(imodel)
    nmems = cmip6models.loc[imodel, "Nmem"]

    for imem in np.arange(0,nmems,1):
        vpdat = cmip6vp_stack.isel(model=imodel, member=imem)
        vpdat_pcnt = cmip6vp_pcnt_stack.isel(model=imodel, member=imem)
        qdat_pcnt = cmip6q_pcnt_stack.isel(model=imodel, member=imem)
        prdat = cmip6pr_stack.isel(model=imodel, member=imem)
        relhumdat = cmip6relhum_stack.isel(model=imodel, member=imem)
        tasdat = cmip6tas_stack.isel(model=imodel, member=imem)

        # Sort into equal area bins according to aridity
        indices = np.argsort(aridity_stack)
        vp_aisort = vpdat[np.array(indices)]
        vp_pcnt_aisort = vpdat_pcnt[np.array(indices)]
        q_pcnt_aisort = qdat_pcnt[np.array(indices)]
        pr_aisort = prdat[np.array(indices)]
        aridity_aisort = aridity_stack[np.array(indices)]
        relhum_aisort = relhumdat[np.array(indices)]
        tas_aisort = tasdat[np.array(indices)]

        # define area weighting for each grid point
        area = np.ones(len(indices))*np.cos(np.deg2rad(aridity_aisort.lat))


        # Split into 50 equal area bins according to aridity index
        cum_area = area.cumsum()/area.sum()
        idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_ai, endpoint=False)[1:])
        area_chunks = np.split(area, idx)
        vp_chunks = np.split(vp_aisort, idx)
        vp_pcnt_chunks = np.split(vp_pcnt_aisort, idx)
        q_pcnt_chunks = np.split(q_pcnt_aisort, idx)
        pr_chunks = np.split(pr_aisort, idx)
        aridity_chunks = np.split(aridity_aisort, idx)
        relhum_chunks = np.split(relhum_aisort, idx)
        tas_chunks = np.split(tas_aisort, idx)

        # Split into 25 equal area bins according to precipitation trend
        for i in np.arange(0,nsplits_ai,1):
            area_bin = area_chunks[i]
            vp_bin = vp_chunks[i]
            vp_pcnt_bin = vp_pcnt_chunks[i]
            q_pcnt_bin = q_pcnt_chunks[i]
            pr_bin = pr_chunks[i]
            aridity_bin = aridity_chunks[i]
            relhum_bin = relhum_chunks[i]
            tas_bin = tas_chunks[i]
            
            indices = np.argsort(pr_bin)
            area_prsort = area_bin[np.array(indices)]
            vp_prsort = vp_bin[np.array(indices)]
            vp_pcnt_prsort = vp_pcnt_bin[np.array(indices)]
            q_pcnt_prsort = q_pcnt_bin[np.array(indices)]
            pr_prsort = pr_bin[np.array(indices)]
            aridity_prsort = aridity_bin[np.array(indices)]
            relhum_prsort = relhum_bin[np.array(indices)]
            tas_prsort = tas_bin[np.array(indices)]

            cum_area = area_prsort.cumsum()/area_prsort.sum()
            idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_pr, endpoint=False)[1:])
            area_prchunks = np.split(area_prsort, idx)
            vp_prchunks = np.split(vp_prsort, idx)
            vp_pcnt_prchunks = np.split(vp_pcnt_prsort, idx)
            q_pcnt_prchunks = np.split(q_pcnt_prsort, idx)
            pr_prchunks = np.split(pr_prsort, idx)
            aridity_prchunks = np.split(aridity_prsort, idx)
            relhum_prchunks = np.split(relhum_prsort, idx)
            tas_prchunks = np.split(tas_prsort, idx)

            for j in np.arange(0,nsplits_pr,1):
                area_cmip6[imodel, imem, i, j] = area_prchunks[j].sum('z')
#                vp_cmip6[imodel, imem, i, j] = vp_prchunks[j].mean('z')
#                vp_cmip6_pcnt[imodel, imem,i,j] = vp_pcnt_prchunks[j].mean('z')
#                q_cmip6_pcnt[imodel, imem, i, j] = q_pcnt_prchunks[j].mean('z') 
#                pr_cmip6[imodel, imem, i, j] = pr_prchunks[j].mean('z')
#                aridity_cmip6[imodel, imem, i, j] = aridity_prchunks[j].mean('z')
#                tas_cmip6[imodel, imem, i, j] = tas_prchunks[j].mean('z')
#                relhum_cmip6[imodel, imem, i, j] = relhum_prchunks[j].mean('z')
                vp_cmip6[imodel, imem, i, j] = meanw(vp_prchunks[j])
                vp_cmip6_pcnt[imodel, imem,i,j] = meanw(vp_pcnt_prchunks[j])
                q_cmip6_pcnt[imodel, imem, i, j] = meanw(q_pcnt_prchunks[j])
                pr_cmip6[imodel, imem, i, j] = meanw(pr_prchunks[j])
                aridity_cmip6[imodel, imem, i, j] = meanw(aridity_prchunks[j])
                tas_cmip6[imodel, imem, i, j] = meanw(tas_prchunks[j])
                relhum_cmip6[imodel, imem, i, j] = meanw(relhum_prchunks[j])






 

dat = xr.merge([area_cmip6, vp_cmip6, vp_cmip6_pcnt, q_cmip6_pcnt, pr_cmip6, aridity_cmip6, relhum_cmip6, tas_cmip6])
dat.to_netcdf("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/fig3/cmip6binned_ai"+
    str(nsplits_ai)+"_pr"+str(nsplits_pr)+"_weighted.nc")







 






