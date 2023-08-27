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
 

stations_t = open("/project/mojave/observations/HadISDH/stations/stationinfo/PosthomogIDPHAt_anoms9120_goodsHadISDH.4.4.0.2021f.txt","r")

latstr=[12,20]
lonstr=[21,30]
statstr=[0,11]

lon_t=[]
lat_t=[]
station_t=[]
for line in stations_t:
    station_t.append(line[statstr[0]:statstr[1]])
    lon_t.append(float(line[lonstr[0]:lonstr[1]]))
    lat_t.append(float(line[latstr[0]:latstr[1]]))   
   
lon_t = xr.DataArray(lon_t, dims=['station_t'], coords=[station_t], name='lon_t')
lat_t = xr.DataArray(lat_t, dims=['station_t'], coords=[station_t], name='lat_t')


 
stations_td = open("/project/mojave/observations/HadISDH/stations/stationinfo/PosthomogPHADPDtd_anoms9120_goodsHadISDH.4.4.0.2021f.txt","r")

latstr=[12,20]
lonstr=[21,30]
statstr=[0,11]

lon_td=[]
lat_td=[]
station_td=[]
for line in stations_td:
    station_td.append(line[statstr[0]:statstr[1]])
    lon_td.append(float(line[lonstr[0]:lonstr[1]]))
    lat_td.append(float(line[latstr[0]:latstr[1]]))

lon_td = xr.DataArray(lon_td, dims=['station_td'], coords=[station_td], name='lon_td')
lat_td = xr.DataArray(lat_td, dims=['station_td'], coords=[station_td], name='lat_td')

   
#-----------------------------------------------------

station_e_set = set(station_e)
station_t_set = set(station_t)
station_td_set = set(station_td)

commonstations = list(station_e_set.intersection(station_t_set, station_td_set))

#----Vapor pressure
basepath="/project/mojave/observations/HadISDH/stations/"
dat=[]
lon=[]
lat=[]
for istation in np.arange(0,len(commonstations),1):
    data_e = xr.open_dataset(basepath+'EDIR/'+commonstations[istation]+"_anoms9120_homog.nc").e_abs
    data_e = data_e.sel(time=slice("1980-01-01","2020-12-31"))
    #dat_e.append(data)

    data_t = xr.open_dataset(basepath+'TDIR/'+commonstations[istation]+"_anoms9120_homog.nc").t_abs
    data_t = data_t.sel(time=slice("1980-01-01","2020-12-31"))
    #dat_e.append(data)

    data_td = xr.open_dataset(basepath+'TDDIR/'+commonstations[istation]+"_anoms9120_homog.nc").td_abs
    data_td = data_td.sel(time=slice("1980-01-01","2020-12-31"))
    #dat_e.append(data) 

    dat.append(xr.merge([data_e, data_t, data_td])) 

    lon.append( ( (lon_e[ lon_e.station_e == commonstations[istation] ]).values)[0])
    lat.append( ( (lat_e[ lon_e.station_e == commonstations[istation] ]).values)[0]) 


data = xr.concat(dat, dim='station')
data.assign_coords({'station':commonstations})
lon = xr.DataArray(lon, dims=['station'], coords=[commonstations], name='lon')
lat = xr.DataArray(lat, dims=['station'], coords=[commonstations], name='lat')

data = xr.merge([data, lon, lat])

data.to_netcdf(pathout+'vp_t_td_HadISDH_stations.nc')

