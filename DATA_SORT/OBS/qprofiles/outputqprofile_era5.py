import importlib
import sys
import xarray as xr
import numpy as np
import xesmf as xe

from qtrendutils import readdata_utils as read
from qtrendutils import calendar_utils as cal
from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

import metpy.interpolate
from math import nan

import dask
dask.config.set(**{'array.slicing.split_large_chunks': True})

sigmaout = np.arange(0.6,1.005,0.005)

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


pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/qprofiles/"


# Read in PS
ps = xr.open_mfdataset("/project/haggis/ERA5/mon/PS/*.nc")
ps = ps.ps
ps = ps.sel(time=slice("1980-01","2020-12"))

# Read in Q
Q = xr.open_mfdataset("/project/haggis/ERA5/mon/Q_mlev/*.nc")
Q = Q.Q
Q = Q.sel(time=slice("1980-01","2020-12"))

# DATA for interpolation grid
cesmdat = xr.open_dataset(
    '../../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

regridder = xe.Regridder(ps, grid_out, 'bilinear', periodic=True, reuse_weights=False, filename="/project/cas/islas/temp/wgtfile.nc")
ps_rg = regridder(ps)

regridder = xe.Regridder(Q, grid_out, 'bilinear', periodic=True, reuse_weights=False, filename="/project/cas/islas/temp/wgtfile.nc")
Q_rg = regridder(Q)


# Read in the hybrid coefficients for ERA5
minlev = np.min(Q.lev).values.astype('int') ; maxlev = np.max(Q.lev).values.astype('int')
erahycoefs = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/ERA5_hycoefs.nc")

hyam = erahycoefs.hyam
hybm = erahycoefs.hybm
hyam = hyam.sel(lev=slice(minlev,maxlev))
hybm = hybm.sel(lev=slice(minlev,maxlev))

states = ['California','Nevada','Utah','Arizona','New Mexico','Colorado']
shpfile="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/shp/gadm36_USA_1.shp"
masksouthwest = shp.maskgen(shpfile, Q_rg,
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

quse = Q_rg.sel(lon=slice(minlon,maxlon), lat=slice(minlat,maxlat))
psuse = ps_rg.sel(lon=slice(minlon,maxlon), lat=slice(minlat,maxlat))
qinterp = lev2sigma(quse,psuse,hyam,hybm,sigmaout,axisinterp="lev", P0=1)

q_sw_temp = avg.cosweightlonlat(qinterp*masksouthwest, 0, 360, -90, 90)
q_sw_temp = q_sw_temp.rename('q_sw')

q_sw_temp.to_netcdf(pathout+'Q_sw_era5.nc')





