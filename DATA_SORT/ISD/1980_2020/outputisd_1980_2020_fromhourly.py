import pandas as pd
import sys
import csv
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan as nan
from qtrendutils import isd_utils as isd
from datetime import datetime
from os.path import exists

ystart=1980 ; yend=2020 ; nyears=yend-ystart+1
timehourly = pd.date_range(str(ystart)+"-01-01",str(yend)+"-12-31", freq="H")
timedaily = pd.date_range(str(ystart)+"-01-01", str(yend)+"-12-31", freq="D")
daystr = xr.DataArray(timehourly.strftime('%Y-%m-%d'), coords=[timehourly], dims=['time'], name='daystr')

datpath="/project/mojave/observations/ISD/hourly/"
statfile="/project/mojave/observations/ISD/global-summary-of-the-day/isd-history.txt"
inventfile="/project/mojave/observations/ISD/global-summary-of-the-day/isd-inventory.csv"
fileout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/ISD/1980_2020/T2M_DEW_ISD_global_dailyfromhourly_QCd_1980_2020.nc"

#---set up the character locations of the required columns fro mthe idf-histor.txt file
datestart=[82,90] 
dateend=[91,99]
latstr=[57,64]
lonstr=[65,73]
usafstr=[0,6]
wbanstr=[7,12]

#---Open the isd-history file and gab out the relevant station information
f = open(statfile,"r")
# skip header
for i in range(22):
    f.readline()

usaf=[] ; wban=[] ; dates=[] ; datee=[] ; lon=[] ; lat=[]
count=0
for line in f:
    dates.append(line[datestart[0]:datestart[1]])
    datee.append(line[dateend[0]:dateend[1]])
    lon.append(line[lonstr[0]:lonstr[1]])
    lat.append(line[latstr[0]:latstr[1]])
    usaf.append(line[usafstr[0]:usafstr[1]])
    wban.append(line[wbanstr[0]:wbanstr[1]])
f.close()

dictstat=[{'wban': wban, 'usaf': usaf, 'lat': lat, 'lon':lon, 'dates':dates, 'datee':datee} 
     for wban, usaf, lat, lon, dates, datee in zip (wban, usaf, lat, lon, dates, datee)]

statinfo = pd.DataFrame.from_dict(dictstat)

print(len(statinfo))

inventory = pd.read_csv(inventfile)

#-----get rid of stations that have white space for lons and lats
#-----or don't cover the full period.  Allowing for 10 days mising at the start and end
bad_ones =  [ (str.isspace(statinfo['lat'][i])) | 
              (str.isspace(statinfo['lon'][i])) | 
              (int(statinfo['dates'][i]) > int(ystart*1e4+10)) | 
              (int(statinfo['datee'][i]) < int(yend*1e4-10)) for i in np.arange(0,len(statinfo),1) ]
bad_ones = pd.DataFrame(bad_ones, index=statinfo.index)
statinfo = statinfo[~bad_ones.values]

statname = statinfo['usaf']+statinfo['wban']
statname = statname.values

print(len(statinfo))



#-----now check if the stations have data for each year.
indices = statinfo.index
for istat in np.arange(0,len(statinfo),1):
    print(istat)
    test = [exists(datpath+str(iyear)+'/'+statname[istat]+'.csv') for iyear in np.arange(ystart,yend+1,1)]
    if (any(i is False for i in test)):
        statinfo = statinfo.drop(indices[istat])

statname = statinfo['usaf']+statinfo['wban']
statname = statname.values

#----Read in each stations data
#----Run the quality controls

usecols=['SOURCE', 'REPORT_TYPE', 'CALL_SIGN', 'QUALITY_CONTROL', 'DATE', 'DEW', 'TMP']
has_data = np.ones(len(statname)).astype(bool)

dat_stations=[]
dat_stations_hourly=[]

for istat in range(0,len(statname),1):
    print(str(istat)+':'+statname[istat])
    

    alldat=[]
    for iyear in range(ystart,yend+1,1):
        fname = datpath+str(iyear)+'/'+statname[istat]+'.csv'

        try:
            data = pd.read_csv(fname, usecols=usecols, low_memory=False)
            data = isd.remove_bad_rows(data)
            data = data[['DATE','TMP','DEW']]
            for varname in ('TMP','DEW'):
                data = isd.remove_bad_vals(data, varname)
            alldat.append(data)
        except:
            has_data[istat] = False
            print("Doesn't have data")
            continue 

    # Check there's data
    if len(alldat) > 0:
        alldat = pd.concat(alldat, ignore_index=True).reset_index()
        alldat = alldat[['DATE','TMP','DEW']]
    else:
        has_data[istat]=False
        continue

    # Give it a time axis
    dt = [datetime.strptime(d, '%Y-%m-%dT%H:%M:%S') for d in alldat['DATE']]
    alldat.index = dt

    #-----check the first and last days are in the record.  If not - pad it with a NaN
    #-----check the first and last days are in the record.  If not - pad it with a NaN
    if ( (alldat.index.year[0] != ystart) | 
         (alldat.index.month[0] != 1) | 
         (alldat.index.day[0] != 1) ):
        day1 = datetime.strptime(str(ystart)+"-01-01", '%Y-%m-%d')
        datday0 = pd.DataFrame(data=[[nan, nan]], index=[day1], columns=['TMP', 'DEW'])
        alldat= pd.concat([datday0, alldat])
    
    if ( (alldat.index.year[len(alldat.index)-1] != yend) | 
         (alldat.index.month[len(alldat.index)-1] != 12) | 
         (alldat.index.day[len(alldat.index)-1] != 31) ):
        dayend = datetime.strptime(str(yend)+"-12-31",'%Y-%m-%d')
        datdayend = pd.DataFrame(data=[[nan, nan]], index=[dayend], columns=['TMP', 'DEW'])
        alldat = pd.concat([alldat, datdayend])


    #-----test for duplicates and keep the first one
    test = alldat.index.duplicated(keep='first')
    alldat = alldat[~test]

    #-----Count the number of observations per day
    resampler = alldat.resample('D')
    nperday = resampler.count()

    #---Do the interpolation
    dathourly=[]
    datdaily=[]

    for varname in ('TMP','DEW'):
        dat = alldat[varname].dropna()
        if (len(dat) > 1):
            datinterp = np.interp(timehourly, dat.index, dat)
            datinterp_xr = xr.DataArray(datinterp, coords=[timehourly], dims='time', name=varname)
            dathourly.append(datinterp_xr)

            datdailyt = datinterp_xr.groupby(daystr).mean('time').rename({'daystr':'time'})
    #       !!!!No dropping < 18 per day right now
   #         datdailyt = datdailyt.where( nperday[varname].values > 18, nan)
            datdailyt['time'] = timedaily

            datdaily.append(datdailyt)
        else:
            has_data[istat] = False

    dathourly = xr.merge(dathourly)
    datdaily = xr.merge(datdaily)

    t2m = datdaily.TMP.rename('t2m')
    dewp = datdaily.DEW.rename('dewp')

    hourly_t2m = dathourly.TMP.rename('t2m')
    hourly_dewp = dathourly.DEW.rename('dewp')
    
    nperday_t2m = xr.DataArray(nperday.TMP.values, coords=[t2m.time], dims=['time'], name='nperday_t2m')
    nperday_dewp = xr.DataArray(nperday.DEW.values, coords=[dewp.time], dims=['time'], name='nperday_dew')

    datout = xr.merge([t2m,dewp, nperday_t2m, nperday_dewp])
    hourly_datout = xr.merge([hourly_t2m, hourly_dewp])

    #---drop Feb 29th
    datout = datout.where( ((datout.time.dt.month != 2) | (datout.time.dt.day != 29)), drop=True)

    dat_stations.append(datout)
    dat_stations_hourly.append(hourly_datout)

dat_stations = xr.concat(dat_stations, dim='station')
dat_stations = dat_stations.assign_coords(station=statname)
lons = xr.DataArray(np.array(statinfo['lon'][has_data].values).astype(float), 
             coords=[dat_stations.station], dims=['station'], name='lons')
lats = xr.DataArray(np.array(statinfo['lat'][has_data].values).astype(float), 
             coords=[dat_stations.station], dims=['station'], name='lats')
dat_stations = xr.merge([dat_stations,lons,lats])
dat_stations.to_netcdf(fileout)



dat_stations_hourly = xr.concat(dat_stations_hourly, dim='station')
dat_stations_hourly = dat_stations_hourly.assign_coords(station=statname)
dat_stations_hourly.to_netcdf("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/ISD/1980_2020/T2M_DEW_ISD_global_hourly_QCd_1980_2020.nc")
