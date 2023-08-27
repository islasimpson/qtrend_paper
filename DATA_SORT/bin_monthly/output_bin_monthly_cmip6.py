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

cmip6models=pd.read_csv('../CMIP6/cmip6csvinfo.csv')

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/bin_monthly/"

#----Land fraction
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
landfrac = landfrac.landfrac
landfrac = landfrac.where(landfrac > 0, nan)
landfrac = landfrac.where( landfrac.lat > -60, nan)


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

aridity_am = aridity_am*landfrac
aridity_monthly = aridity_monthly*landfrac

# mask the monthly aridity indices to have limite the aridity index to 10 to prevent
# blow ups when dividing by small numbers
aridity_monthly = aridity_monthly.where(
      ( (aridity_monthly < 10) & (np.isfinite(aridity_monthly)) & ~np.isnan(aridity_monthly)), 10)



#----GPCC precipitation trend (used to remove grid points that don't have GPCC data)
#gpcc_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/'+
#    'DATA_SORT/1980_2020_trends/prtrends_GPCC.nc')
#gpcc_map['lon'] = landfrac.lon ; gpcc_map['lat'] = landfrac.lat
gpcc_clim = xr.open_dataset(
      "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/"
       +"gpccclim_monthly.nc")
gpcc_clim = gpcc_clim.mean('month', skipna=False)
gpcc_clim['lon'] = landfrac.lon ; gpcc_clim['lat'] = landfrac.lat
gpcc_clim_stack = gpcc_clim.stack(z=('lon','lat'))


aridity_stack = aridity_am.stack(z=('lon','lat'))

idrop = np.argwhere( np.isnan(np.array(aridity_stack[:])) |
                     ~np.isfinite(np.array(aridity_stack[:])) |
                     np.isnan(np.array(gpcc_clim_stack.pr) ) )[:,0]


#gpcc_trend_stack = gpcc_map.stack(z=('lon','lat'))
#idrop = np.argwhere( np.isnan(np.array(aridity_stack[:])) |
#                     ~np.isfinite(np.array(aridity_stack[:])) |
#                     np.isnan(np.array(gpcc_trend_stack.pr) ) )[:,0]


#-----CMIP6 data
cmip6vp_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/vptrends_CMIP6_monthly.nc')
cmip6vp_map['lon'] = landfrac.lon ; cmip6vp_map['lat'] = landfrac.lat
cmip6vp_map = cmip6vp_map*landfrac

cmip6pr_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/prtrends_CMIP6_monthly.nc')*86400.
cmip6pr_map['lon'] = landfrac.lon ; cmip6pr_map['lat'] = landfrac.lat
cmip6pr_map = cmip6pr_map*landfrac

cmip6q_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/qtrends_CMIP6_monthly.nc')
cmip6q_map['lon'] = landfrac.lon ; cmip6q_map['lat'] = landfrac.lat
cmip6q_map = cmip6q_map*landfrac

cmip6tas_map = xr.open_dataset('/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_trends/tastrends_CMIP6_monthly.nc')
cmip6tas_map['lon'] = landfrac.lon ; cmip6tas_map['lat'] = landfrac.lat
cmip6tas_map = cmip6tas_map*landfrac

# CMIP6 climatologies (for % calculation)
#cmip6clim = xr.open_dataset(
#  "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/"+
#  "vpclims_CMIP6_monthly.nc")
#cmip6clim['lon'] = landfrac.lon ; cmip6clim['lat'] = landfrac.lat
#cmip6clim_vp_1990 = cmip6clim.vp*landfrac
#cmip6clim_q_1990 = cmip6clim.q*landfrac
#cmip6clim_tas_1990 = cmip6clim.tas*landfrac

cmip6clim_vp_1990 = xr.open_dataset(
   '/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/'+
   'vpclims_CMIP6_monthly.nc').vp
cmip6clim_vp_1990['lon'] = landfrac.lon ; cmip6clim_vp_1990['lat'] = landfrac.lat
cmip6clim_vp_1990 = cmip6clim_vp_1990*landfrac

cmip6clim_q_1990 = xr.open_dataset(
   '/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/'+
   'qclims_CMIP6_monthly.nc').q
cmip6clim_q_1990['lon'] = landfrac.lon ; cmip6clim_q_1990['lat'] = landfrac.lat
cmip6clim_q_1990 = cmip6clim_q_1990*landfrac

cmip6clim_tas_1990 = xr.open_dataset(
   '/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_1990_clims/'+
   'tasclims_CMIP6_monthly.nc').tas
cmip6clim_tas_1990['lon'] = landfrac.lon ; cmip6clim_tas_1990['lat'] = landfrac.lat
cmip6clim_tas_1990 = cmip6clim_tas_1990*landfrac

cmip6vp_map_pcnt = (cmip6vp_map / cmip6clim_vp_1990)*100.
cmip6q_map_pcnt = (cmip6q_map / cmip6clim_q_1990)*100.

#cmip6clim_2020 = xr.open_dataset(
#  "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/"+
#  "vpclims_CMIP6_monthly.nc")
#cmip6clim_2020['lon'] = landfrac.lon ; cmip6clim_2020['lat'] = landfrac.lat
#cmip6clim_vp_2020 = cmip6clim.vp*landfrac
#cmip6clim_q_2020 = cmip6clim.q*landfrac
#cmip6clim_tas_2020 = cmip6clim.tas*landfrac

cmip6clim_vp_2020 = xr.open_dataset(
   '/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/'+
   'vpclims_CMIP6_monthly.nc').vp
cmip6clim_vp_2020['lon'] = landfrac.lon ; cmip6clim_vp_2020['lat'] = landfrac.lat
cmip6clim_vp_2020 = cmip6clim_vp_2020*landfrac

cmip6clim_q_2020 = xr.open_dataset(
   '/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/'+
   'qclims_CMIP6_monthly.nc').q
cmip6clim_q_2020['lon'] = landfrac.lon ; cmip6clim_q_2020['lat'] = landfrac.lat
cmip6clim_q_2020 = cmip6clim_q_2020*landfrac

cmip6clim_tas_2020 = xr.open_dataset(
   '/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/'+
   'tasclims_CMIP6_monthly.nc').tas
cmip6clim_tas_2020['lon'] = landfrac.lon ; cmip6clim_tas_2020['lat'] = landfrac.lat
cmip6clim_tas_2020 = cmip6clim_tas_2020*landfrac

cmip6clim_pr_2020 = xr.open_dataset(
  "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/1980_2020_clims/"+
  "prclims_CMIP6_monthly.nc").pr
cmip6clim_pr_2020['lon'] = landfrac.lon ; cmip6clim_pr_2020['lat'] = landfrac.lat
cmip6clim_pr_2020 = cmip6clim_pr_2020*86400.



#cmip6vp_map_stack = cmip6vp_map.stack(z=('member','model'))
#cmip6vp_map_stack = cmip6vp_map_stack.dropna('z')
#
#cmip6pr_map_stack = cmip6pr_map.stack(z=('member','model'))
#cmip6pr_map_stack = cmip6pr_map_stack.dropna('z')

#----CMIP6 regression coefficients
regcoefs = xr.open_dataset(pathout+'cmip6_monthly_regcoefs.nc')
vppredict = regcoefs.aglobal + regcoefs.bglobal*cmip6pr_map.pr
vppredict = vppredict.transpose("model","member","month","lat","lon")
trenddif = cmip6vp_map.vp - np.array(vppredict)


#----Now doing the binning
aridity_monthly_stack = aridity_monthly.stack(loc=('lon','lat')).rename('aridity_monthly')
aridity_am_stack = aridity_am.stack(loc=('lon','lat')).rename('aridity_am')
vppredict_stack = vppredict.stack(loc=('lon','lat')).rename({'vppredict'})
trenddif_stack = trenddif.stack(loc=('lon','lat')).rename('vptrenddif')
vp_stack = cmip6vp_map.stack(loc=('lon','lat'))
q_stack = cmip6q_map.stack(loc=('lon','lat'))
pr_stack = cmip6pr_map.stack(loc=('lon','lat'))
vp_pcnt_stack = cmip6vp_map_pcnt.stack(loc=('lon','lat'))
q_pcnt_stack = cmip6q_map_pcnt.stack(loc=('lon','lat'))
tas_stack = cmip6tas_map.stack(loc=('lon','lat'))
vp_clim_1990_stack = cmip6clim_vp_1990.stack(loc=('lon','lat'))
q_clim_1990_stack = cmip6clim_q_1990.stack(loc=('lon','lat'))
tas_clim_1990_stack = cmip6clim_tas_1990.stack(loc=('lon','lat'))
vp_clim_2020_stack = cmip6clim_vp_2020.stack(loc=('lon','lat'))
q_clim_2020_stack = cmip6clim_q_2020.stack(loc=('lon','lat'))
tas_clim_2020_stack = cmip6clim_tas_2020.stack(loc=('lon','lat'))
pr_clim_2020_stack = cmip6clim_pr_2020.stack(loc=('lon','lat'))


aridity_monthly_stack = aridity_monthly_stack.drop_isel(loc=idrop)
aridity_am_stack = aridity_am_stack.drop_isel(loc=idrop)
vppredict_stack = vppredict_stack.drop_isel(loc=idrop)
trenddif_stack = trenddif_stack.drop_isel(loc=idrop)
vp_stack = vp_stack.drop_isel(loc=idrop)
q_stack = q_stack.drop_isel(loc=idrop)
pr_stack = pr_stack.drop_isel(loc=idrop)
vp_pcnt_stack = vp_pcnt_stack.drop_isel(loc=idrop)
q_pcnt_stack = q_pcnt_stack.drop_isel(loc=idrop)
tas_stack = tas_stack.drop_isel(loc=idrop)
vp_clim_1990_stack = vp_clim_1990_stack.drop_isel(loc=idrop)
q_clim_1990_stack = q_clim_1990_stack.drop_isel(loc=idrop)
tas_clim_1990_stack = tas_clim_1990_stack.drop_isel(loc=idrop)
vp_clim_2020_stack = vp_clim_2020_stack.drop_isel(loc=idrop)
q_clim_2020_stack = q_clim_2020_stack.drop_isel(loc=idrop)
tas_clim_2020_stack = tas_clim_2020_stack.drop_isel(loc=idrop)
pr_clim_2020_stack = pr_clim_2020_stack.drop_isel(loc=idrop)


# test for monthly aridity indices that are still nan's because pet=0, set them to the maximum finite ai
aridity_monthly_stack = aridity_monthly_stack.where( np.isfinite(aridity_monthly_stack), nan)
maxai = np.nanmax(aridity_monthly_stack)
aridity_monthly_stack = aridity_monthly_stack.where(
           ~np.isnan(aridity_monthly_stack) & np.isfinite(aridity_monthly_stack), maxai)


lons = vp_stack.lon.rename('lons')
lats = vp_stack.lat.rename('lats')

nsplits_ai = 30 # 30 aridity bins

vp_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vp')

q_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='q')

pr_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='pr')

vppredict_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vppredict')
 
trenddif_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='trenddif')

aridity_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='aridity')

vp_pcnt_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vp_pcnt')

q_pcnt_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='q_pcnt')

tas_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='tas')

tas_clim_1990_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='tas_clim_1990')

vp_clim_1990_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vp_clim_1990')

q_clim_1990_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='q_clim_1990')

tas_clim_2020_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='tas_clim_2020')

vp_clim_2020_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='vp_clim_2020')

q_clim_2020_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='q_clim_2020')

pr_clim_2020_cmip6 = xr.DataArray(
    np.zeros([cmip6vp_map.model.size, cmip6vp_map.member.size, nsplits_ai, 12])*nan,
    coords=[cmip6vp_map.model, cmip6vp_map.member, np.arange(0,nsplits_ai,1), np.arange(0,12,1)],
    dims=['model','member','ai','month'], name='pr_clim_2020')




lons = lons.reset_index('loc')
lats = lats.reset_index('loc')


#----cosine weighted mean over data points
def meanw(dat, lats):
    meanw = np.sum(dat*np.cos(np.deg2rad(np.array(lats)))) / np.sum(np.cos(np.deg2rad(np.array(lats))))
    return meanw


lons = lons.reset_coords( ('lon','lat'), drop=True)
lats = lats.reset_coords( ('lon','lat'), drop=True)



#----loop over models and members
for imodel in np.arange(0,cmip6vp_map.model.size,1):
    print(imodel)
    nmems = cmip6models.loc[imodel, "Nmem"]

    for imem in np.arange(0,nmems,1):
        vpdat = vp_stack.vp.isel(model=imodel, member=imem)
        qdat = q_stack.q.isel(model=imodel, member=imem)
        prdat = pr_stack.pr.isel(model=imodel, member=imem)
        vppredictdat = vppredict_stack.isel(model=imodel, member=imem)
        trenddifdat = trenddif_stack.isel(model=imodel, member=imem)
        vp_pcnt_dat = vp_pcnt_stack.vp.isel(model=imodel, member=imem)
        q_pcnt_dat = q_pcnt_stack.q.isel(model=imodel, member=imem)
        tas_dat = tas_stack.tas.isel(model=imodel, member=imem)
        vp_clim_1990_dat = vp_clim_1990_stack.isel(model=imodel, member=imem)
        q_clim_1990_dat = q_clim_1990_stack.isel(model=imodel, member=imem)
        tas_clim_1990_dat = tas_clim_1990_stack.isel(model=imodel, member=imem)
        vp_clim_2020_dat = vp_clim_2020_stack.isel(model=imodel, member=imem)
        q_clim_2020_dat = q_clim_2020_stack.isel(model=imodel, member=imem)
        tas_clim_2020_dat = tas_clim_2020_stack.isel(model=imodel, member=imem)
        pr_clim_2020_dat = pr_clim_2020_stack.isel(model=imodel, member=imem)


        #-----Sort along the month axis according to aridity
        aridity_sorted = np.take_along_axis(aridity_monthly_stack.values, 
              aridity_monthly_stack.argsort(axis=0), axis=0)
        vp_sorted = np.take_along_axis(vpdat.values, aridity_monthly_stack.argsort(axis=0),
                     axis=0)
        q_sorted = np.take_along_axis(qdat.values, aridity_monthly_stack.argsort(axis=0),
                     axis=0)
        vppredict_sorted = np.take_along_axis(vppredictdat.values, aridity_monthly_stack.argsort(axis=0),
                     axis=0)
        pr_sorted = np.take_along_axis(prdat.values, aridity_monthly_stack.argsort(axis=0),
                     axis=0)
        trenddif_sorted = np.take_along_axis(trenddifdat.values, aridity_monthly_stack.argsort(axis=0),
                     axis=0)
        vp_pcnt_sorted = np.take_along_axis(vp_pcnt_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        q_pcnt_sorted = np.take_along_axis(q_pcnt_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        tas_sorted = np.take_along_axis(tas_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        vp_clim_1990_sorted = np.take_along_axis(vp_clim_1990_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        q_clim_1990_sorted = np.take_along_axis(q_clim_1990_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        tas_clim_1990_sorted = np.take_along_axis(tas_clim_1990_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        vp_clim_2020_sorted = np.take_along_axis(vp_clim_2020_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        q_clim_2020_sorted = np.take_along_axis(q_clim_2020_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        tas_clim_2020_sorted = np.take_along_axis(tas_clim_2020_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)
        pr_clim_2020_sorted = np.take_along_axis(pr_clim_2020_dat.values, aridity_monthly_stack.argsort(axis=0),
                            axis=0)




        aridity_sorted = xr.DataArray(aridity_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='aridity')
        vp_sorted = xr.DataArray(vp_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='vp')
        q_sorted = xr.DataArray(q_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='q')
        vppredict_sorted = xr.DataArray(vppredict_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='vppredict')
        pr_sorted = xr.DataArray(pr_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='pr')
        trenddif_sorted = xr.DataArray(trenddif_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='trenddif')
        vp_pcnt_sorted = xr.DataArray(trenddif_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='vp_pcnt')
        q_pcnt_sorted = xr.DataArray(q_pcnt_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='q_pcnt')
        tas_sorted = xr.DataArray(tas_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='tas')
        vp_clim_1990_sorted = xr.DataArray(vp_clim_1990_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='vp_clim_1990')
        q_clim_1990_sorted = xr.DataArray(q_clim_1990_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='q_clim_1990')
        tas_clim_1990_sorted = xr.DataArray(tas_clim_1990_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='tas_clim_1990')
        vp_clim_2020_sorted = xr.DataArray(vp_clim_2020_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='vp_clim_2020')
        q_clim_2020_sorted = xr.DataArray(q_clim_2020_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='q_clim_2020')
        tas_clim_2020_sorted = xr.DataArray(tas_clim_2020_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='tas_clim_2020')
        pr_clim_2020_sorted = xr.DataArray(pr_clim_2020_sorted, coords=aridity_monthly_stack.coords,
              dims=aridity_monthly_stack.dims, name='pr_clim_2020')




 
        alldat = xr.merge([aridity_sorted, vp_sorted, q_sorted, vppredict_sorted, pr_sorted,
                           trenddif_sorted, aridity_am_stack, lons, lats,
                           vp_pcnt_sorted, q_pcnt_sorted, tas_sorted, vp_clim_1990_sorted, 
                           q_clim_1990_sorted,tas_clim_1990_sorted, vp_clim_2020_sorted,
                           q_clim_2020_sorted, tas_clim_2020_sorted, pr_clim_2020_sorted])

        alldatsort = alldat.sortby(alldat.aridity_am)

        area = np.ones(alldatsort.lats.size)*np.cos(np.deg2rad(alldatsort.lats))

        # split into equal area bins according to aridity index
        cum_area = area.cumsum() / area.sum()
        idx = np.searchsorted(cum_area, np.linspace(0,1,nsplits_ai, endpoint=False)[1:])

        lons_chunks = np.split(alldatsort.lons, idx)
        lats_chunks = np.split(alldatsort.lats, idx)
        vp_chunks = np.split(alldatsort.vp, idx, axis=1)
        q_chunks = np.split(alldatsort.q, idx, axis=1)
        vppredict_chunks = np.split(alldatsort.vppredict, idx, axis=1)
        trenddif_chunks = np.split(alldatsort.trenddif, idx, axis=1)
        pr_chunks = np.split(alldatsort.pr, idx, axis=1)
        aridity_chunks = np.split(alldatsort.aridity, idx, axis=1)
        vp_pcnt_chunks = np.split(alldatsort.vp_pcnt, idx, axis=1)
        q_pcnt_chunks = np.split(alldatsort.q_pcnt, idx, axis=1)
        tas_chunks = np.split(alldatsort.tas, idx, axis=1)
        vp_clim_1990_chunks = np.split(alldatsort.vp_clim_1990, idx, axis=1)
        q_clim_1990_chunks = np.split(alldatsort.q_clim_1990, idx, axis=1)
        tas_clim_1990_chunks = np.split(alldatsort.tas_clim_1990, idx, axis=1)
        vp_clim_2020_chunks = np.split(alldatsort.vp_clim_2020, idx, axis=1)
        q_clim_2020_chunks = np.split(alldatsort.q_clim_2020, idx, axis=1)
        tas_clim_2020_chunks = np.split(alldatsort.tas_clim_2020, idx, axis=1)
        pr_clim_2020_chunks = np.split(alldatsort.pr_clim_2020, idx, axis=1)





        for i in np.arange(0,nsplits_ai, 1):
            for j in np.arange(0,12,1):
                vp_cmip6[imodel,imem,i,j] = meanw( np.array(vp_chunks[i])[j,:], np.array(lats_chunks[i]))
                q_cmip6[imodel,imem,i,j] = meanw( np.array(q_chunks[i])[j,:], np.array(lats_chunks[i]))
                vppredict_cmip6[imodel,imem,i,j] = meanw( np.array(vppredict_chunks[i])[j,:], np.array(lats_chunks[i]))
                trenddif_cmip6[imodel,imem,i,j] = meanw( np.array(trenddif_chunks[i])[j,:], np.array(lats_chunks[i]))
                pr_cmip6[imodel,imem,i,j] = meanw( np.array(pr_chunks[i])[j,:], np.array(lats_chunks[i]))
                aridity_cmip6[imodel,imem,i,j] = meanw( np.array(aridity_chunks[i])[j,:], np.array(lats_chunks[i]))
                vp_pcnt_cmip6[imodel,imem,i,j] = meanw( np.array(vp_pcnt_chunks[i])[j,:], np.array(lats_chunks[i]))
                q_pcnt_cmip6[imodel,imem,i,j] = meanw( np.array(q_pcnt_chunks[i])[j,:], np.array(lats_chunks[i]))
                tas_cmip6[imodel, imem, i, j] = meanw( np.array(tas_chunks[i])[j,:], np.array(lats_chunks[i]))
                vp_clim_1990_cmip6[imodel,imem,i,j] = meanw( np.array(vp_clim_1990_chunks[i])[j,:], np.array(lats_chunks[i]))
                q_clim_1990_cmip6[imodel,imem,i,j] = meanw( np.array(q_clim_1990_chunks[i])[j,:], np.array(lats_chunks[i]))
                tas_clim_1990_cmip6[imodel,imem,i,j] = meanw( np.array(tas_clim_1990_chunks[i])[j,:], np.array(lats_chunks[i]))
                vp_clim_2020_cmip6[imodel,imem,i,j] = meanw( np.array(vp_clim_2020_chunks[i])[j,:], np.array(lats_chunks[i]))
                q_clim_2020_cmip6[imodel,imem,i,j] = meanw( np.array(q_clim_2020_chunks[i])[j,:], np.array(lats_chunks[i]))
                tas_clim_2020_cmip6[imodel,imem,i,j] = meanw( np.array(tas_clim_2020_chunks[i])[j,:], np.array(lats_chunks[i]))
                pr_clim_2020_cmip6[imodel,imem,i,j] = meanw( np.array(pr_clim_2020_chunks[i])[j,:], np.array(lats_chunks[i]))




 
datout = xr.merge([vp_cmip6, q_cmip6, vppredict_cmip6, trenddif_cmip6, pr_cmip6, aridity_cmip6,
           vp_pcnt_cmip6, q_pcnt_cmip6, vp_clim_1990_cmip6, q_clim_1990_cmip6, tas_clim_1990_cmip6, 
           vp_clim_2020_cmip6, q_clim_2020_cmip6, tas_clim_2020_cmip6, pr_clim_2020_cmip6,tas_cmip6])
datout.to_netcdf(pathout+'CMIP6_binned_monthly.nc')                
