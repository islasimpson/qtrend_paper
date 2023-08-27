import importlib
import xarray as xr
import numpy as np
import sys
import pandas as pd
import warnings
import math
import xesmf as xe
from glob import glob

import matplotlib.pyplot as plt

from qtrendutils import readdata_utils as read

warnings.filterwarnings('ignore')

histpath="/project/cmip5/historical/Amon/"
rcppath="/project/cmip5/rcp85/Amon/"

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP5/pr/"

cmip5models=pd.read_csv('cmip5csvinfo.csv')

ybegp = 1950 ; monbegp = 1 ; yendp = 2014 ; monendp = 12 # dates for Past period
ybegf = 2015 ; monbegf = 1 ; yendf = 2021 ; monendf = 12 # dates for Future period

# total number of months (used for checking)
nmonthsp = (yendp-ybegp-1)*12 + (12-monbegp+1) + monendp
nmonthsf = (yendf-ybegf-1)*12 + (12-monbegf+1) + monendf
nyearsf = yendf - ybegf + 1
nyearsp = yendp - ybegp + 1

# set up date names
datebegp=str(ybegp)+"-"+str(monbegp).zfill(2)
dateendp=str(yendp)+"-"+str(monendp).zfill(2)
datebegf=str(ybegf)+"-"+str(monbegf).zfill(2)
dateendf=str(yendf)+"-"+str(monendf).zfill(2)

# DATA for interpolation grid
cesmdat = xr.open_dataset(
    '../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

models = cmip5models['Model']
nmodels = models.size

wgtfile=pathout+"wgtfile.nc"

for index, modname in models.iteritems():
    print(modname)
    nmems = cmip5models.loc[index, "Nmem"]

    wgtreuse=False

    for imem in range(1,nmems+1,1):
        memstr="r"+str(imem)+"i1p1"
        # check if a special order is needed
        changeorder = cmip5models.loc[index,"specialorder"]
        if (type(changeorder) == str):
            changeordernp = np.array(changeorder.split(","))
            memstr = "r"+str(changeordernp[imem-1])+"i1p1"

        histdir_pr = glob(histpath+"/pr/"+modname+"/"+memstr)   
        histdir_pr = histdir_pr[0]
    
        rcpdir_pr = glob(rcppath+"/pr/"+modname+"/"+memstr)
        rcpdir_pr = rcpdir_pr[0]
    
        pr_hist = read.read_sfc(histdir_pr+"/*.nc", datebegp, dateendp)
        pr_rcp = read.read_sfc(rcpdir_pr+"/*.nc", datebegf, dateendf)
    
        pr = xr.concat([pr_hist, pr_rcp], dim="time", join="override")
    
        regridder = xe.Regridder(pr, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                        filename=wgtfile)
    
        wgtreuse = True
        pr_rg = regridder(pr.pr)
    
        if (imem == 1):
            pr_out = xr.DataArray(np.zeros([nmems, pr.time.size, pr_rg.lat.size, pr_rg.lon.size]),
                      coords = [np.arange(1,nmems+1,1), pr.time, pr_rg.lat, pr_rg.lon],
                      dims=['member','time','lat','lon'], name='pr')
    
        pr_out[imem-1,:,:,:] = pr_rg
    
    pr_out.to_netcdf(pathout+"pr_"+modname+".nc")

