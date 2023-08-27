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
            if ( ~np.isnan(landfrac[ilat,ilon])):
                a_t, b_t = linfit.linfit_xy(cmip6pr_map_stack.pr_cmip6.isel(month=imon, lon=ilon, lat=ilat), cmip6vp_map_stack.vp_cmip6.isel(month=imon, lon=ilon, lat=ilat))
                a_global[imon,ilat,ilon] = a_t ; b_global[imon,ilat,ilon] = b_t

a_global = xr.DataArray(a_global, coords=[cmip6pr_map_stack.month, cmip6pr_map_stack.lat, cmip6pr_map_stack.lon], dims=['month','lat','lon'], name='aglobal')
b_global = xr.DataArray(b_global, coords=[cmip6pr_map_stack.month, cmip6pr_map_stack.lat, cmip6pr_map_stack.lon], dims=['month','lat','lon'], name='bglobal')

a_global.to_netcdf(pathout+'cmip6_monthly_regcoefs.nc')
b_global.to_netcdf(pathout+'cmip6_monthly_regcoefs.nc', mode='a')
