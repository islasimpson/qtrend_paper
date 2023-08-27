import importlib
import sys
import xarray as xr
import numpy as np

from qtrendutils import lensread_utils as lens
from qtrendutils import readdata_utils as read
from qtrendutils import calendar_utils as cal
from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

import metpy.interpolate
from math import nan

memnames = lens.lens2memnamegen(100)

sigmaout=np.arange(0.4,1.005,0.005)

def lev2sigma(dat, ps, hyam, hybm, sigmaout, axisinterp="lev", P0=1e5):
    """ Interpolate from hybrid sigma coordinates to sigma coordinates"""


    pre = hyam*P0 + hybm*ps
    sigma = pre/ps

    dat = dat.transpose("lev",...)
    sigma = sigma.transpose("lev",...)

    def verticalinterp(dat, sigma):
        return metpy.interpolate.interpolate_1d(sigmaout, sigma, dat, axis=0)

    output = xr.apply_ufunc(
             verticalinterp,
             dat,
             sigma,
             input_core_dims=[[axisinterp],[axisinterp]],
             output_core_dims=[['sigma']],
             vectorize=True,
             dask = "parallelized",
             output_dtypes=[dat.dtype],
             dask_gufunc_kwargs = {"output_sizes": {"sigma":len(sigmaout)}})

    output['sigma'] = sigmaout

    return(output)

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CESM2/LENS2/Qprofiles/"
qdir="/project/mojave/cesm2/LENS/atm/month_1/Q/"
psdir="/project/mojave/cesm2/LENS/atm/month_1/PS/"

for imem in range(0,len(memnames),1):
    print(memnames[imem])
    qfiles = qdir+"*-"+memnames[imem]+"*.nc"
    q = xr.open_mfdataset(qfiles)
    q = read.fixcesmtime(q)
    q = q.sel(time=slice("1980-01","2022-12"))
    hyam = q.hyam
    hybm = q.hybm
    q = q.Q

    psfiles = psdir+"*-"+memnames[imem]+"*.nc"
    ps = xr.open_mfdataset(psfiles)
    ps = read.fixcesmtime(ps)
    ps = ps.sel(time=slice("1980-01","2022-12"))
    ps = ps.PS

    if (imem == 0):
        states = ['California','Nevada','Utah','Arizona','New Mexico','Colorado']
        dat = xr.open_dataset(
         '/project/mojave/cesm2/LENS/atm/month_1/Q/'+\
         'b.e21.BSSP370smbb.f09_g17.LE2-1301.020.cam.h0.Q.209501-210012.nc')
        shpfile="/project/cas/islas/shapefiles/usa/gadm36_USA_1.shp"
        masksouthwest = shp.maskgen(shpfile, dat,
                 ['California','Nevada','Utah','Arizona','New Mexico','Colorado'])

        latarray = np.tile(masksouthwest.lat, masksouthwest.lon.size)
        latarray = np.reshape(latarray, [masksouthwest.lon.size, masksouthwest.lat.size])
        latarray = np.moveaxis(latarray, 1, 0)

        lonarray = np.tile(masksouthwest.lon, masksouthwest.lat.size)
        lonarray = np.reshape(lonarray, [masksouthwest.lat.size, masksouthwest.lon.size])

        lonsw = np.where( masksouthwest == 1, lonarray,  nan)
        latsw = np.where( masksouthwest == 1, latarray, nan)

        minlon = np.nanmin(lonsw) ; maxlon = np.nanmax(lonsw)
        minlat = np.nanmin(latsw) ; maxlat = np.nanmax(latsw)

        masksouthwest = masksouthwest.sel(lon=slice(minlon,maxlon), lat=slice(minlat,maxlat))

    quse = q.sel(lon=slice(minlon,maxlon), lat=slice(minlat, maxlat))
    psuse = ps.sel(lon=slice(minlon,maxlon), lat=slice(minlat, maxlat))
    qinterp = lev2sigma(quse, psuse, hyam, hybm, sigmaout, axisinterp="lev")

    q_sw_temp = avg.cosweightlonlat(qinterp*masksouthwest, 0, 360, -90, 90)

    q_sw_temp = q_sw_temp.rename('q_sw')

    q_sw_temp.to_netcdf(pathout+'Q_sw_lens2_'+memnames[imem]+'.nc')
