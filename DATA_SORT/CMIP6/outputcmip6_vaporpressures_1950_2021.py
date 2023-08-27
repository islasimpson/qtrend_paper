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
from qtrendutils import humiditycalcs as qcalcs

warnings.filterwarnings('ignore')

histpath="/project/mojave/cmip6/historical/Amon/"
ssppath="/project/mojave/cmip6/ssp585/Amon/"

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/CMIP6/vapor_pressures/"

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

        histdir_huss = glob(histpath+"/huss/"+modname+"/"+memstr+"/*")
        histdir_huss = histdir_huss[0]

        sspdir_huss = glob(ssppath+"/huss/"+modname+"/"+memstr+"/*")
        sspdir_huss = sspdir_huss[0]

        histdir_tas = glob(histpath+"/tas/"+modname+"/"+memstr+"/*")
        histdir_tas = histdir_tas[0]

        sspdir_tas = glob(ssppath+"/tas/"+modname+"/"+memstr+"/*")
        sspdir_tas = sspdir_tas[0]

        histdir_ps = glob(histpath+"/ps/"+modname+"/"+memstr+"/*")
        #!!!!Dealing with IITM-ESM for now
        if (modname == 'IITM-ESM' ):
            histdir_ps = glob('/project/cas/islas/CMIP6/historical/Amon/ps/IITM-ESM/')
        histdir_ps = histdir_ps[0]

        sspdir_ps = glob(ssppath+"/ps/"+modname+"/"+memstr+"/*")
        sspdir_ps = sspdir_ps[0]

        huss_hist = read.read_sfc(histdir_huss+'/*.nc', datebegp, dateendp)
        tas_hist = read.read_sfc(histdir_tas+'/*.nc', datebegp, dateendp)
        ps_hist = read.read_sfc(histdir_ps+'/*.nc', datebegp, dateendp)

        huss_ssp = read.read_sfc(sspdir_huss+'/*.nc', datebegf, dateendf)
        tas_ssp = read.read_sfc(sspdir_tas+'/*.nc', datebegf, dateendf)
        ps_ssp = read.read_sfc(sspdir_ps+'/*.nc', datebegf, dateendf)

        huss = xr.concat([huss_hist, huss_ssp], dim="time", join="override")
        tas = xr.concat([tas_hist, tas_ssp], dim="time", join="override")
        ps = xr.concat([ps_hist, ps_ssp], dim="time", join="override")

        regridder = xe.Regridder(ps, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                  filename=wgtfile)

        wgtreuse = True
        huss_rg = regridder(huss.huss)
        tas_rg = regridder(tas.tas)
        ps_rg = regridder(ps.ps)

        #--Saturation vapor pressure
        svp = qcalcs.calcsvp(tas_rg)
   
        #--Vapor pressure
        vp = qcalcs.calcvpfromhuss(huss_rg, ps_rg)

        #--Vapor pressure deficit
        vpd = svp - vp

        #--Saturation specific humidity
        sq = qcalcs.calcsq(tas_rg, ps_rg)

        #--Relative humidity
        qrel = (huss_rg / sq)*100.
        qrel = qrel.assign_attrs({'units':'%'})


        if (imem == 1):
            svp_out = xr.DataArray(np.zeros([nmems, svp.time.size, ps_rg.lat.size, ps_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), svp.time, ps_rg.lat, ps_rg.lon],
                  dims=['member','time','lat','lon'], name='svp',
                  attrs={'long_name':'Saturation Vapor Pressure', 'units':'hPa'})
            vp_out = xr.DataArray(np.zeros([nmems, vp.time.size, ps_rg.lat.size, ps_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), vp.time, ps_rg.lat, ps_rg.lon],
                  dims=['member','time','lat','lon'], name='vp',
                  attrs={'long_name':'Vapor Pressure', 'units':'hPa'})
            vpd_out = xr.DataArray(np.zeros([nmems, vpd.time.size, ps_rg.lat.size, ps_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), vpd.time, ps_rg.lat, ps_rg.lon],
                  dims=['member','time','lat','lon'], name='vpd',
                  attrs={'long_name':'Vapor Pressure Deficit', 'units':'hPa'})
            qrel_out = xr.DataArray(np.zeros([nmems, qrel.time.size, ps_rg.lat.size, ps_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), qrel.time, ps_rg.lat, ps_rg.lon],
                  dims=['member','time','lat','lon'], name='relhum',
                  attrs={'long_name':'Relative Humidity', 'units':'%'})
            tas_out = xr.DataArray(np.zeros([nmems, tas_rg.time.size, tas_rg.lat.size, tas_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), tas_rg.time, tas_rg.lat, tas_rg.lon],
                  dims=['member','time','lat','lon'], name='tas',
                  attrs={'long_name':'Near surface air temperature', 'units':'K'})
            q_out = xr.DataArray(np.zeros([nmems, huss_rg.time.size, huss_rg.lat.size, huss_rg.lon.size]),
                  coords = [np.arange(1,nmems+1,1), huss_rg.time, huss_rg.lat, huss_rg.lon],
                  dims=['member','time','lat','lon'], name='q',
                  attrs={'long_name':'Near surface specific humidity', 'units':'kg/kg'})




        svp_out[imem-1,:,:,:] = svp
        vp_out[imem-1,:,:,:] = vp
        vpd_out[imem-1,:,:,:] = vpd
        qrel_out[imem-1,:,:,:] = qrel
        tas_out[imem-1,:,:,:] = tas_rg
        q_out[imem-1,:,:,:] = huss_rg

    svp_out.to_netcdf(pathout+"vaporpressures_"+modname+".nc")
    vp_out.to_netcdf(pathout+"vaporpressures_"+modname+".nc", mode="a")
    vpd_out.to_netcdf(pathout+"vaporpressures_"+modname+".nc", mode="a")
    qrel_out.to_netcdf(pathout+"vaporpressures_"+modname+".nc", mode="a")
    tas_out.to_netcdf(pathout+"vaporpressures_"+modname+".nc", mode="a")
    q_out.to_netcdf(pathout+"vaporpressures_"+modname+".nc", mode="a")
