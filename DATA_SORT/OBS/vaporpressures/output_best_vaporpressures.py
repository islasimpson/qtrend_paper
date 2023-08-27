import xarray as xr
import numpy as np
import xesmf as xe

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs
from qtrendutils import calendar_utils as cal

pathout = "/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"
pathT = "//project/mojave/observations/OBS-TAS/"

#datT = read.read_sfc(pathT+'best.tas.185001-202112.nc', "1980-01","2020-12")
datT = xr.open_dataset(pathT+'best.tas.185001-202112.nc')
datT = datT.tas
datT = datT+273.15


time = cal.YYYYMM2date(datT.time)
datT['time'] = time
datT = datT.sel(time=slice("1980-01-01","2020-12-31"))

#-----Regridding T and TD before the vp and svp calculation for consistency with CMIP
wgtfile=pathout+'wgtfile.nc'
cesmdat = xr.open_dataset(
     '../../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})
wgtreuse=False
regridder = xe.Regridder(datT, grid_out, 'bilinear', periodic=True, reuse_weights=wgtreuse,
                         filename=wgtfile)

datT_rg = regridder(datT)

#--Saturation vapor pressure
svp = qcalcs.calcsvp(datT_rg)
svp = svp.rename('svp')

t2m = datT_rg
t2m = t2m.rename('t2m')
t2m = t2m.assign_attrs({'units':'K', 'long_name':'K'})

svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})

svp.to_netcdf(pathout+'vaporpressures_BEST.nc')
t2m = t2m.to_netcdf(pathout+'vaporpressures_BEST.nc', mode='a')

