import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/HadISDH/"

latstr=[12,20]
lonstr=[21,30]
statstr=[0,11]

#---Read in the stations for the individual variables
stations_e = open("/project/mojave/observations/HadISDH/stations/stationinfo/PosthomogIDPHAe_anoms9120_goodsHadISDH.4.4.0.2021f.txt","r")

latstr=[12,20]
lonstr=[21,30]
statstr=[0,11]

lon_e=[]
lat_e=[]
station_e=[]
for line in stations_e:
    station_e.append(line[statstr[0]:statstr[1]])
    lon_e.append(float(line[lonstr[0]:lonstr[1]]))
    lat_e.append(float(line[latstr[0]:latstr[1]]))

lon_e = xr.DataArray(lon_e, dims=['station_e'], coords=[station_e], name='lon_e')
lat_e = xr.DataArray(lat_e, dims=['station_e'], coords=[station_e], name='lat_e')
 

#-----------------------------------------------------


#----Vapor pressure
basepath="/project/mojave/observations/HadISDH/stations/"
dat=[]
lon=[]
lat=[]
for istation in np.arange(0,len(station_e),1):
    data_e = xr.open_dataset(basepath+'EDIR/'+station_e[istation]+"_anoms9120_homog.nc").e_abs
    data_e = data_e.sel(time=slice("1980-01-01","2020-12-31"))
    #dat_e.append(data)

    dat.append(data_e) 

    lon.append((lon_e[istation]).values)
    lat.append((lat_e[istation]).values) 


data = xr.concat(dat, dim='station')
data = data.rename('vp')
data.assign_coords({'station':station_e})
lon = xr.DataArray(lon, dims=['station'], coords=[station_e], name='lon')
lat = xr.DataArray(lat, dims=['station'], coords=[station_e], name='lat')

data = xr.merge([data, lon, lat])


data.to_netcdf(pathout+'vp_direct_HadISDH_stations.nc')

