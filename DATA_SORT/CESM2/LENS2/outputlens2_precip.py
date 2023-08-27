import xarray as xr
import numpy as np
import sys

from qtrendutils import humiditycalcs as qcalcs
from qtrendutils import lensread_utils as lens
from qtrendutils import readdata_utils as read

memnames = lens.lens2memnamegen(100)

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/precip/"
precc_dir="/project/mojave/cesm2/LENS/atm/month_1/PRECC/"
precl_dir="/project/mojave/cesm2/LENS/atm/month_1/PRECL/"

for imem in range(0,len(memnames),1):
    print(memnames[imem])
    preccfiles = precc_dir+"*-"+memnames[imem]+"*.nc"
    preclfiles = precl_dir+"*-"+memnames[imem]+"*.nc"

    precc = xr.open_mfdataset(preccfiles)
    precc = read.fixcesmtime(precc)
    precc = precc.sel(time=slice("1980-01","2020-12"))
    precc = precc.PRECC

    precl = xr.open_mfdataset(preclfiles)
    precl = read.fixcesmtime(precl)
    precl = precl.sel(time=slice("1980-01","2020-12"))
    precl = precl.PRECL

    prect = (precc+precl)*86400.*1000.

    prect = prect.rename('pr')
    prect = prect.assign_attrs({'units':'mm/day', 'long_name':'Total precipitation'})

    prect.to_netcdf(pathout+'pr_LENS2_'+memnames[imem]+'.nc')
