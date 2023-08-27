import xarray as xr
import numpy as np
import sys

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs
from qtrendutils import lensread_utils as lens

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/GOGA2/"
preccdir="/project/mojave/cesm2/f.e21.FHIST_FSSP370_BGC.f09_f09.ersstv5.goga/atm/month_1/PRECC/"
precldir="/project/mojave/cesm2/f.e21.FHIST_FSSP370_BGC.f09_f09.ersstv5.goga/atm/month_1/PRECL/"

for imem in np.arange(0,10,1):
    memstr = str(imem+1).zfill(2)
    print(memstr)

    precc = xr.open_mfdataset(preccdir+'*.goga*'+str(imem+1).zfill(2)+'.cam*.nc')
    precl = xr.open_mfdataset(precldir+'*.goga*'+str(imem+1).zfill(2)+'.cam*.nc')

    precc = read.fixcesmtime(precc)
    precl = read.fixcesmtime(precl)

    prect = precc.PRECC + precl.PRECL

    prect = prect*86400.*1000.
    prect = prect.rename('pr')

    #--Total precipitation
    prect = prect.assign_attrs({'units':'mm/day', 'long_name':'Total precipitation'})

    prect.to_netcdf(pathout+'pr_GOGA2_'+memstr+'.nc')





