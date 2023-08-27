import xarray as xr
import numpy as np
import sys

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs
from qtrendutils import lensread_utils as lens

memnames = lens.lens2memnamegen(100)

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/vaporpressures/"
trefhtdir="/project/mojave/cesm2/LENS/atm/month_1/TREFHT/"
psdir="/project/mojave/cesm2/LENS/atm/month_1/PS/"
qrefhtdir="/project/mojave/cesm2/LENS/atm/month_1/QREFHT/"


for imem in range(0,len(memnames),1):
    print(memnames[imem])

    trefhtfiles = trefhtdir+"*-"+memnames[imem]+"*.nc"
    qrefhtfiles = qrefhtdir+"*-"+memnames[imem]+"*.nc"
    psfiles = psdir+"*-"+memnames[imem]+"*.nc"

    trefht = xr.open_mfdataset(trefhtfiles)
    trefht = read.fixcesmtime(trefht)
    trefht = trefht.sel(time=slice("1980-01","2020-12"))
    trefht = trefht.TREFHT

    qrefht = xr.open_mfdataset(qrefhtfiles)
    qrefht = read.fixcesmtime(qrefht)
    qrefht = qrefht.sel(time=slice("1980-01","2020-12")) 
    qrefht = qrefht.QREFHT

    ps = xr.open_mfdataset(psfiles)
    ps = read.fixcesmtime(ps)
    ps = ps.sel(time=slice("1980-01","2020-12"))
    ps = ps.PS


    #--Saturation vapor pressure
    svp = qcalcs.calcsvp(trefht)
    svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
    svp = svp.rename('svp')

    #--Vapor pressure
    vp = qcalcs.calcvpfromhuss(qrefht, ps)
    vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})
    vp = vp.rename('vp')

    #--Vapor pressure deficit
    vpd = svp - vp
    vpd = vpd.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure Deficit'})
    vpd = vpd.rename('vpd')

    #--Saturation specific humidity
    sq = qcalcs.calcsq(trefht, ps)
    sq = sq.assign_attrs({'units':'hPa','long_name':'Saturation Specific Humidity'})
    sq = sq.rename('svp')


    #--Relative humidity
    qrel = (qrefht / sq)*100.
    qrel = qrel.assign_attrs({'units':'%', 'long_name':'Relative Humidity'})
    qrel = qrel.rename('relhum')    

    t2m = trefht
    t2m = t2m.assign_attrs({'units':'K', 'long_name':'Near surface temperature'})

    qrefht = qrefht.rename('q')


    datout = xr.merge([svp, vp, vpd, sq, qrel, t2m,qrefht], compat='override')

    datout.to_netcdf(pathout+'vaporpressures_LENS2_'+memnames[imem]+'.nc')


