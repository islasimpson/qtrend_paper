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

histpath="/project/mojave/cmip6/historical/Amon/"
ssppath="/project/mojave/cmip6/ssp585/Amon/"

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP6/pr/"

cmip6models=pd.read_csv('cmip6csvinfo.csv')

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


models = cmip6models['Model']
nmodels = models.size

wgtfile=pathout+"wgtfile.nc"

for index, modname in models.iteritems():
    print(modname)
    nmems = cmip6models.loc[index, "Nmem"]
    ftype = cmip6models.loc[index, "ftypep"]
    ptype = cmip6models.loc[index, "ptype"]

    wgtreuse=False

    for imem in range(1, nmems+1,1):
        if (math.isnan(ptype)):
            memstr="r"+str(imem)+"i1p1f"+str(int(ftype))
        else:
            memstr="r"+str(imem)+"i1p"+str(int(ptype))+"f"+str(int(ftype))

        # check if s apecial order is needed
        changeorder = cmip6models.loc[index, "specialorder"]
        if (type(changeorder) == str):
            changeordernp = np.array(changeorder.split(","))
            if (math.isnan(ptype)):
                memstr = "r"+str(changeordernp[imem-1])+"i1p1f"+str(int(ftype))
            else:
                memstr = "r"+str(changeordernp[imem-1])+"i1p"+str(int(ptype))+"f"+str(int(ftype))

        histdir_pr = glob(histpath+"/pr/"+modname+"/"+memstr+"/*")
        histdir_pr = histdir_pr[0]

        sspdir_pr = glob(ssppath+"/pr/"+modname+"/"+memstr+"/*")
        sspdir_pr = sspdir_pr[0]

        pr_hist = read.read_sfc(histdir_pr+"/*.nc", datebegp, dateendp)
        pr_ssp = read.read_sfc(sspdir_pr+"/*.nc", datebegf, dateendf)

        pr = xr.concat([pr_hist, pr_ssp], dim="time", join="override")

        regridder = xe.Regridder(pr, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                  filename=wgtfile)

        wgtreuse = True
        pr_rg = regridder(pr.pr)

        if (modname == 'CIESM'):
            pr_rg = pr_rg*1000.

        if (imem == 1):
            pr_out = xr.DataArray(np.zeros([nmems, pr.time.size, pr_rg.lat.size, pr_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), pr.time, pr_rg.lat, pr_rg.lon],
                  dims=['member','time','lat','lon'], name='pr')

        pr_out[imem-1,:,:,:] = pr_rg

    pr_out.to_netcdf(pathout+"pr_"+modname+".nc")




