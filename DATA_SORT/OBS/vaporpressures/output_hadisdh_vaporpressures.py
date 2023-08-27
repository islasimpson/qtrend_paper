import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from math import nan

from qtrendutils import humiditycalcs as qcalcs

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"

hadisdh = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/HadISDH/"+
   "vp_t_td_HadISDH_stations.nc")

badvp = xr.where( np.isnan(hadisdh.e_abs), 1, 0)
badtd = xr.where( np.isnan(hadisdh.td_abs), 1, 0)
badt = xr.where( np.isnan(hadisdh.t_abs), 1, 0)
#
#nbadvp = badvp.sum("time")
#nbadtd = badtd.sum("time")
#nbadt = badt.sum("time")
#
#nbadvp = nbadvp.rename('nbadvp')
#nbadtd = nbadtd.rename('nbadtd')
#nbadt = nbadt.rename('nbadt')

#----drop the lons and lats and add them back in later
lon = hadisdh.lon
lat = hadisdh.lat
isd = hadisdh.drop(['lon','lat'])

#----Read in ERA5 to obtain the surface pressure at the closest grid point
era5ps = xr.open_mfdataset("/project/mojave/observations/ERA5/PS/*.nc")

time1=str(hadisdh.time[0].dt.year.values)+'-'+str(hadisdh.time[0].dt.month.values).zfill(2)
time2=str(hadisdh.time[hadisdh.time.size-1].dt.year.values)+'-'+ \
      str(hadisdh.time[hadisdh.time.size - 1].dt.month.values)

era5ps = era5ps.sel(time=slice(time1,time2))


#----Flip the longitudes of HadISDH to go from 0 to 360
newlons = xr.where(hadisdh['lon'] < 0, hadisdh['lon'] + 360, hadisdh['lon'])
hadisdh['lon'] = newlons

#----Find the ERA5 ps at the equivalent grid point to the stations
ps = xr.DataArray(np.zeros([hadisdh.station.size, hadisdh.time.size]),
       coords=[hadisdh.station, hadisdh.time], 
       dims=['station','time'], name='ps')
for istation in np.arange(0,hadisdh.station.size,1):
    print(istation)
    lonval = hadisdh.lon.isel(station=istation).values
    latval = hadisdh.lat.isel(station=istation).values
    lonargmin = np.argmin( np.abs( era5ps.lon.values -  lonval))
    psuse = era5ps.isel(lon=lonargmin)
    latargmin = np.argmin( np.abs( era5ps.lat.values - latval))
    psuse = psuse.isel(lat=latargmin)
    ps[istation,:] = psuse.ps.values

hadisdh = xr.merge([hadisdh, ps])

#-----Get dew point T and T in K
t2dew = hadisdh.td_abs + 273.15

t2m = hadisdh.t_abs + 273.15

ps = hadisdh.ps

#-----Saturation vapor pressure
svp = qcalcs.calcsvp(t2m)
svp = svp.rename('svp')

#-----Vapor pressure
vp = qcalcs.calcsvp(t2dew)
vp = vp.rename('vp')

vp_fromdataset = hadisdh.e_abs
vp_fromdataset = vp_fromdataset.rename('vp_fromhadisdh')


#-----Vapor pressure deficit
vpd = svp - vp
vpd = vpd.rename('vpd')

#-----Saturation specific humidity
sq = qcalcs.calcsq(t2m, ps)

#-----Specific humidity
q = qcalcs.calcsq(t2dew,ps)

#-----Relative humidity
qrel = (q / sq)*100.
qrel = qrel.rename('relhum')

qrel = qrel.assign_attrs({'units':'%', 'long_name':'Relative humidity'})
svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation vapor pressure'})
vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})
vpd = vpd.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure Deficit'})

lon = hadisdh.lon
lat = hadisdh.lat

qrel.to_netcdf(pathout+'vaporpressures_HadISDH.nc')
svp.to_netcdf(pathout+'vaporpressures_HadISDH.nc', mode='a')
vp.to_netcdf(pathout+'vaporpressures_HadISDH.nc', mode='a')
vpd.to_netcdf(pathout+'vaporpressures_HadISDH.nc', mode='a')
vp_fromdataset.to_netcdf(pathout+'vaporpressures_HadISDH.nc', mode='a')
lon.to_netcdf(pathout+'vaporpressures_HadISDH.nc', mode='a')
lat.to_netcdf(pathout+'vaporpressures_HadISDH.nc', mode='a')




