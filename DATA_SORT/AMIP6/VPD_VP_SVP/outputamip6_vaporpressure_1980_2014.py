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

histpath="/project/mojave/cmip6/amip/Amon/"

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/AMIP6/vapor_pressures/"

cmip6models=pd.read_csv('/home/islas/python/qtrend_paper/DATA_SORT/AMIP6/amip6csvinfo.csv')

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

        histdir_tas = glob(histpath+"/tas/"+modname+"/"+memstr+"/*")
        histdir_tas = histdir_tas[0]

        histdir_ps = glob(histpath+"/ps/"+modname+"/"+memstr+"/*")
        histdir_ps = histdir_ps[0]

#        if ( (modname == 'CAS-ESM2-0') & (memstr == 'r1i1p1f1')):
        if ( modname == 'CAS-ESM2-0'):
            huss = read.read_sfc(histdir_huss+"/*.nc", datebegp, dateendp)
            tas = read.read_sfc(histdir_tas+"/*_1979*.nc", datebegp, dateendp)
            ps = read.read_sfc(histdir_ps+"/*.nc", datebegp, dateendp)
        elif (modname == 'TaiESM1'):
#        if (modname == 'TaiESM1'):
            huss = read.read_sfc(histdir_huss+"/*.nc", datebegp, dateendp)
            tas = read.read_sfc(histdir_tas+"/*_197901-201412.nc", datebegp, dateendp)
            ps = read.read_sfc(histdir_ps+"/*_197901-201412.nc", datebegp, dateendp)
        else:
            huss = read.read_sfc(histdir_huss+"/*.nc", datebegp, dateendp)
            tas = read.read_sfc(histdir_tas+"/*.nc", datebegp, dateendp)
            ps = read.read_sfc(histdir_ps+"/*.nc", datebegp, dateendp)
            
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
