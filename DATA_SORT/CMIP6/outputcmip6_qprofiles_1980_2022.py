import importlib
import sys
import xarray as xr
import numpy as np
import pandas as pd
import math
from glob import glob

from qtrendutils import readdata_utils as read
from qtrendutils import calendar_utils as cal
from qtrendutils import shapefile_utils as shp
from qtrendutils import averaging_utils as avg

import metpy.interpolate
from math import nan

import xesmf as xe
import warnings
warnings.filterwarnings('ignore')

sigmaout=np.arange(0.4,1.005,0.005)

def lev2sigma(dat, ps, sigmaout, axisinterp="plev"):
    """ Interpolate from hybrid sigma coordinates to sigma coordinates"""


    sigma = dat.plev/ps

    dat = dat.transpose("plev",...)
    sigma = sigma.transpose("plev",...)

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


histpath="/project/mojave/cmip6/historical/Amon/"
ssppath="/project/mojave/cmip6/ssp585/Amon/"

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP6/Qprofiles/"

cmip6models = pd.read_csv('cmip6csvinfo_hus.csv')

ybegp = 1980 ; monbegp = 1 ; yendp = 2014 ; monendp = 12 # dates for Past period
ybegf = 2015 ; monbegf = 1 ; yendf = 2022 ;  monendf = 12 # dates for Future period

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

        print(memstr)
        histdir_hus = glob(histpath+"/hus/"+modname+"/"+memstr+"/*")
        histdir_hus = histdir_hus[0]

        sspdir_hus = glob(ssppath+"/hus/"+modname+"/"+memstr+"/*")
        sspdir_hus = sspdir_hus[0]

        histdir_ps = glob(histpath+"/ps/"+modname+"/"+memstr+"/*")


        if (modname == 'IITM-ESM' ):
            histdir_ps = glob('/project/cas/islas/CMIP6/historical/Amon/ps/IITM-ESM/')
        histdir_ps = histdir_ps[0]

        sspdir_ps = glob(ssppath+"/ps/"+modname+"/"+memstr+"/*")
        sspdir_ps = sspdir_ps[0]


        hus_hist = read.read_sfc(histdir_hus+'/*.nc', datebegp, dateendp)
        ps_hist = read.read_sfc(histdir_ps+'/*.nc', datebegp, dateendp)

        hus_ssp = read.read_sfc(sspdir_hus+'/*.nc', datebegf, dateendf)
        ps_ssp = read.read_sfc(sspdir_ps+'/*.nc', datebegf, dateendf)

        hus = xr.concat([hus_hist, hus_ssp], dim="time", join="override")
        ps = xr.concat([ps_hist, ps_ssp], dim="time", join="override")

        # setting the fill values to NaN's.  
        hus = hus.where( (np.abs(hus.hus) < 1e10), nan)

        #----Fill in the missing data
        hus = hus.hus
        levindices = xr.DataArray(np.arange(0,hus.plev.size,1), coords=[hus.plev], 
                                  dims=['plev'], name='levindices')
        levindices4d = hus*0 + levindices
        levbot = xr.where( levindices4d == levindices4d.min('plev'), 1, nan)
        filldat4d = hus*levbot
        filldat = xr.DataArray(np.zeros([hus.time.size, hus.lat.size, hus.lon.size]),
                     coords=[hus.time, hus.lat, hus.lon], dims=['time','lat','lon'], name='filldat')
        for ilev in np.arange(0,hus.plev.size,1):
            datlev = filldat4d.isel(plev=ilev)
            filldat = filldat.where( np.isnan(datlev), datlev)
        filldat = filldat.drop_vars("plev")
        filldat = filldat.expand_dims(dim={"plev":hus.plev.values}, axis=1)
        hus = hus.where( ~np.isnan(hus), filldat)
        #----End fill in missing data

        if (imem == 1):
            wgtreuse = False
        else:
            wgtreuse = True

        regridder_ps = xe.Regridder(ps, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                  filename=wgtfile)
        regridder_hus = xe.Regridder(hus, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                  filename=wgtfile)

       
        q = regridder_hus(hus)
        ps = regridder_ps(ps.ps)

        states = ['California','Nevada','Utah','Arizona','New Mexico','Colorado']
        shpfile="/project/cas/islas/shapefiles/usa/gadm36_USA_1.shp"
        masksouthwest = shp.maskgen(shpfile, ps,
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

        qinterp = lev2sigma(quse, psuse, sigmaout, axisinterp="plev")

        q_sw_temp = avg.cosweightlonlat(qinterp*masksouthwest, 0, 360, -90, 90)
        q_sw_temp = q_sw_temp.rename('q_sw')

        if (imem == 1):
            q_sw_out = xr.DataArray(np.zeros([nmems, q_sw_temp.time.size, len(sigmaout)]),
                       coords=[np.arange(1,nmems+1,1), q_sw_temp.time, sigmaout],
                       dims=['member','time','sigmaout'], name='q_sw')

        q_sw_out[imem-1,:,:] = q_sw_temp[:,:]

    q_sw_out.to_netcdf(pathout+'Q_sw_'+modname+'.nc')
