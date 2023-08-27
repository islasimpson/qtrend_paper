import importlib
import xarray as xr
import numpy as np
import sys
import pandas as pd
import warnings
import math
import xesmf as xe

import matplotlib.pyplot as plt
from CASutils import mapplot_utils as maps
from CASutils import readdata_utils as read
from glob import glob

importlib.reload(maps)
importlib.reload(read)

warnings.filterwarnings('ignore')

histpath="/project/mojave/cmip6/amip/Amon/"
pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/AMIP6/PRECIP/"

cmip6models = pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/AMIP6/amip6csvinfo.csv')

ybegp = 1980 ; monbegp = 1 ; yendp = 2014 ; monendp = 12 # dates for Past period

# total number of months (used for checking)
nmonthsp = (yendp-ybegp-1)*12 + (12-monbegp+1) + monendp
nyearsp = yendp - ybegp + 1

# set up date names
datebegp=str(ybegp)+"-"+str(monbegp).zfill(2)
dateendp=str(yendp)+"-"+str(monendp).zfill(2)

# DATA for interpolation grid
cesmdat = xr.open_dataset(
    '/fs/cgd/csm/inputdata/atm/cam2/inic/fv/cami-mam3_0000-01-01_0.9x1.25_L32_c141031.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

models = cmip6models['Model']
nmodels = models.size

wgtfile = pathout+"wgtfile.nc"

for index, modname in models.iteritems():
    print(modname)
    nmemsp = cmip6models.loc[index, "Nmem"]
    ftype = cmip6models.loc[index, "ftypep"]
    ptype = cmip6models.loc[index, "ptype"]

    wgtreuse=False

    for imem in range(1, nmemsp+1,1):
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

        if ( (modname == 'CAS-ESM2-0') & (memstr == "r1i1p1f1")):
            pr = read.read_sfc(histdir_pr+"/*_1979*.nc", datebegp, dateendp)
        elif (modname == "TaiESM1"):
            pr = read.read_sfc(histdir_pr+"/*_197901-201412.nc", datebegp, dateendp)
        else:
            pr = read.read_sfc(histdir_pr+"/*.nc", datebegp, dateendp)


        regridder = xe.Regridder(pr, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                  filename=wgtfile)

        wgtreuse = True
        pr_rg = regridder(pr.pr)


        if (imem == 1):
            pr_out = xr.DataArray(np.zeros([nmemsp, pr.time.size, pr_rg.lat.size, pr_rg.lon.size]),
                  coords = [np.arange(1,nmemsp+1,1), pr.time, pr_rg.lat, pr_rg.lon],
                  dims=['member','time','lat','lon'], name='pr')

        pr_out[imem-1,:,:,:] = pr_rg*86400.
        pr_out.attrs['units'] = 'mm/day'

    pr_out.to_netcdf(pathout+'precip_'+modname+'.nc')



