import xarray as xr
import numpy as np
import sys

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs
from qtrendutils import lensread_utils as lens

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/GOGA2/vaporpressures_climps/"
trefhtdir="/project/mojave/cesm2/f.e21.FHIST_FSSP370_BGC.f09_f09.ersstv5.goga/atm/month_1/TREFHT/"
psdir="/project/mojave/cesm2/f.e21.FHIST_FSSP370_BGC.f09_f09.ersstv5.goga/atm/month_1/PS/"
qrefhtdir="//project/mojave/cesm2/f.e21.FHIST_FSSP370_BGC.f09_f09.ersstv5.goga/atm/month_1/QREFHT/"

for imem in np.arange(0,10,1):
    memstr = str(imem+1).zfill(2)
    print(memstr)

    trefht = xr.open_mfdataset(trefhtdir+'*.goga*'+str(imem+1).zfill(2)+'.cam*.nc')
    qrefht = xr.open_mfdataset(qrefhtdir+'*.goga*'+str(imem+1).zfill(2)+'.cam*.nc')
    ps = xr.open_mfdataset(psdir+'*.goga*'+str(imem+1).zfill(2)+'.cam*.nc')

    trefht = read.fixcesmtime(trefht)
    qrefht = read.fixcesmtime(qrefht)
    ps = read.fixcesmtime(ps)

    trefht = trefht.sel(time=slice("1980-01-01","2020-12-31"))
    qrefht = qrefht.sel(time=slice("1980-01-01","2020-12-31"))
    ps = ps.sel(time=slice("1980-01-01","2020-12-31"))

    trefht = trefht.TREFHT
    qrefht = qrefht.QREFHT
    ps = ps.PS

    psclim = ps.groupby('time.month').mean('time')
    psclim_expand = psclim.expand_dims(dim={'year':np.arange(0,ps.time.size/12.,1)})
    ps_stack = psclim_expand.stack(z=("year","month"))
    ps_stack = ps_stack.transpose("z","lat","lon")
    ps_stack = ps_stack.rename(z="time")
    ps_stack['time'] = ps.time

    #--Saturation vapor pressure
    svp = qcalcs.calcsvp(trefht)
    svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
    svp = svp.rename('svp')

    #--Vapor pressure
    vp = qcalcs.calcvpfromhuss(qrefht, ps_stack)
    vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})
    vp = vp.rename('vp')

    #--Vapor pressure deficit
    vpd = svp - vp
    vpd = vpd.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure Deficit'})
    vpd = vpd.rename('vpd')

    #--Saturation specific humidity
    sq = qcalcs.calcsq(trefht, ps_stack)
    sq = sq.assign_attrs({'units':'hPa','long_name':'Saturation Specific Humidity'})
    sq = sq.rename('svp')


    #--Relative humidity
    qrel = (qrefht / sq)*100.
    qrel = qrel.assign_attrs({'units':'%', 'long_name':'Relative Humidity'})
    qrel = qrel.rename('relhum')

    datout = xr.merge([svp, vp, vpd, sq, qrel], compat='override')

    datout.to_netcdf(pathout+'vaporpressures_GOGA2_'+memstr+'.nc')





