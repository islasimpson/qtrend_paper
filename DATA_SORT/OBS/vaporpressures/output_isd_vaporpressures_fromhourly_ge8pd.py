#---Output the vapor pressures, VPD and relative humidity for ISD
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from math import nan

from qtrendutils import humiditycalcs as qcalcs

isd = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/ISD/1980_2020/T2M_DEW_ISD_global_dailyfromhourly_QCd_1980_2020.nc", use_cftime=True)

#----Set the days that have less than 4 observations to NaN's
t2m = isd.t2m
dewp = isd.dewp
t2m = t2m.where(isd.nperday_t2m >= 8, nan)
dewp = dewp.where(isd.nperday_dew >= 8, nan)
isd['t2m'] = t2m
isd['dewp'] = dewp


#----drop the lons and lats and add them back in later
lons = isd.lons
lats = isd.lats
isd = isd.drop(['lons','lats'])


#----Check for the number of dew point NaN's in a month
monyearstr = xr.DataArray(isd.indexes['time'].strftime('%Y-%m'), coords=isd.time.coords, name='monyearstr')
dewp = isd.dewp
dewp999s = xr.where( ((dewp > 998) | np.isnan(dewp)),1,0)
nbaddays = dewp999s.groupby(monyearstr).sum("time")
nbaddays = nbaddays.rename('nbaddays')
nbaddays = nbaddays.rename(monyearstr='time')


#----Drop the stations that have more than 48 bad months in the record
#dropstations=[]
#for istation in np.arange(0,nbaddays.station.size,1):
#    test = nbaddays.isel(station=istation)
#    test = test.where(test > 15, drop=True)
#    if (test.size > 48):
#        badstation = nbaddays.station.isel(station=istation).values
#        badstation = np.array_str(badstation)
#        dropstations.append(badstation)
#
#isd = isd.drop_sel(station=dropstations)
#nbaddays = nbaddays.drop_sel(station=dropstations)

    
#---convert all those with a dew point temperature greater than 998 to Nan
isd = isd.where(isd.dewp < 998, nan)
#isd = xr.merge([isd,lon,lat])

#------------------------------------------------------------

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"

#----Calculate the monthly mean
timeout = isd.time.groupby(monyearstr).mean('time')
timeout = timeout.rename(monyearstr='time')

isd_monthly = isd.groupby(monyearstr).mean('time', skipna=True)
isd_monthly = isd_monthly.rename(monyearstr='time')
isd_monthly['time'] = timeout
nbaddays['time'] = timeout

isd_monthly = xr.merge([isd_monthly, lons, lats])
#--------



#----Read in ERA5 to obtain the surface pressure at the closest grid point
era5ps = xr.open_mfdataset("/project/mojave/observations/ERA5/PS/*.nc")
time1 = str(isd_monthly.time.dt.year.isel(time=0).values)+'-'+str(isd_monthly.time.dt.month.isel(time=0).values).zfill(2)
time2 = str(isd_monthly.time.dt.year.isel(time=isd_monthly.time.size-1).values)+'-'+\
        str(isd_monthly.time.dt.month.isel(time=isd_monthly.time.size-1).values).zfill(2)
era5ps = era5ps.sel(time=slice(time1,time2))
#--------

#----Flip the longitudes of IDF to go from 0 to 360
newlons = xr.where(isd_monthly['lons'] < 0, isd_monthly['lons'] + 360., isd_monthly['lons'])
isd_monthly['lons'] = newlons
#--------

#----Find ERA5 ps at the equivalent grid point to the stations
ps = xr.DataArray(np.zeros([isd_monthly.time.size, isd_monthly.station.size]), coords=[isd_monthly.time, isd_monthly.station],
                    dims=['time','station'], name='ps')
for istation in np.arange(0,isd_monthly.station.size,1):
    print(istation)
    lonval = isd_monthly.lons.isel(station=istation).values
    latval = isd_monthly.lats.isel(station=istation).values
    lonargmin = np.argmin( np.abs(era5ps.lon.values - lonval)) 
    psuse = era5ps.isel(lon=lonargmin)
    latargmin = np.argmin( np.abs(era5ps.lat.values - latval))
    psuse = psuse.isel(lat=latargmin)
    ps[:,istation] = psuse.ps.values

isd_monthly = xr.merge([isd_monthly, ps])
#---------

#---------Get dew point T and T in K.  Omit bad data
t2dew = isd_monthly.dewp
#t2dew = t2dew.where( t2dew < 1000, nan)
#t2dew = (t2dew-32.)/1.8 + 273.15
t2dew = t2dew + 273.15

t2m = isd_monthly.t2m
#t2m = (t2m-32)/1.8 + 273.15
t2m = t2m+273.15

ps = isd_monthly.ps

#---------Saturation vapor pressure
svp = qcalcs.calcsvp(t2m)
svp = svp.rename('svp')

#---------Vapor pressure
vp = qcalcs.calcsvp(t2dew)
vp = vp.rename('vp')

#---------Vapor pressure deficit
vpd = svp - vp
vpd = vpd.rename('vpd')

#---------Saturation specific humidity
sq = qcalcs.calcsq(t2m,ps)

#---------Specific humidity
q = qcalcs.calcsq(t2dew, ps)

#---------Relative humidity
qrel = (q / sq)*100.
qrel = qrel.rename('relhum')

qrel = qrel.assign_attrs({'units':'%', 'long_name':'Relative Humidity'})
svp = svp.assign_attrs({'units':'hPa', 'long_name':'Saturation Vapor Pressure'})
vp = vp.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure'})
vpd = vpd.assign_attrs({'units':'hPa', 'long_name':'Vapor Pressure Deficit'})

lon = isd_monthly.lons
lat = isd_monthly.lats

qrel.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc')
svp.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc', mode='a')
vp.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc', mode='a')
vpd.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc', mode='a')
lon.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc', mode='a')
lat.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc', mode='a')
nbaddays.to_netcdf(pathout+'vaporpressures_ISD_fromhourly_ge8pd.nc', mode='a')




