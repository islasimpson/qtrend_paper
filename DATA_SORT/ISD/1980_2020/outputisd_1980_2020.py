import pandas as pd
import sys 
import csv
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from CASutils import calendar_utils as cal
from CASutils import mapplot_utils as maps
from math import nan as nan

ystart=1980 ; yend=2020 ; nyears=yend-ystart+1

### Setting up the output calendar dates
timeout = pd.date_range(start=str(ystart)+"-01-01",end=str(yend)+"-12-31")
### Remove Feb 29th
timeout = timeout[~((timeout.month == 2) & (timeout.day == 29))]

monstrings=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
datpath="/project/mojave/observations/ISD/global-summary-of-the-day/archive/"
statfile="/project/mojave/observations/ISD/global-summary-of-the-day/isd-history.txt"
inventfile="/project/mojave/observations/ISD/global-summary-of-the-day/isd-inventory.csv"
fileout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/ISD/1980_2020/T2M_DEW_ISD_global_1980_2020.nc"

#---Set up the character locations of the required columns from the idf-history.txt file
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

#sys.exit()
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

inventory = pd.read_csv(inventfile)

statuse=[]
statname=[]
usaf=[]
wban=[]
latstation=[]
lonstation=[]

for key in dictstat:
    try:
        latflt=float(key['lat'])
    except:
        latflt=-9999.
        
    try: 
        lonflt=float(key['lon'])
    except:
        lonflt=-9999.
        
    try:
        usafval = int(key['usaf'])
    except:
        usafval = key['usaf']
        
    inventdat = inventory.loc[inventory['USAF']==usafval]
    statyears = inventdat.loc[(inventdat['YEAR'] >= ystart) & (inventdat['YEAR'] <= yend)]
    statyears = statyears['YEAR']
    
    datebegflt = float(key['dates'])
    dateendflt = float(key['datee'])
    
    if (  (datebegflt < ystart*10000) and (dateendflt >= (yend*10000+1)) and (len(statyears) == nyears)):
        statuse.append(key)
        statname.append(key['usaf']+key['wban'])
        usaf.append(key['usaf'])
        wban.append(key['wban'])
        latstation.append(latflt)
        lonstation.append(lonflt)


datevals=pd.date_range(start="1979-01-01",end="1979-12-31")
m = np.array(datevals.month)
mm = np.char.zfill(m.astype(str),2)
d = np.array(datevals.day)
dd = np.char.zfill(d.astype(str),2)
datestrings=[mmm + "-" +  ddd  for  mmm, ddd in zip(mm, dd) ]

# Read in the data for stations that will be used
t2m = np.empty([len(statname), nyears*365])
t2m[:,:] = nan
dewp = np.empty([len(statname), nyears*365])
dewp[:,:] = nan
stp = np.empty([len(statname), nyears*365])
stp[:,:] = nan

for istat in range(0,len(statname),1):
    print(str(istat)+':'+statname[istat])

    usafval = int(usaf[istat])
    wbanval = int(wban[istat])
    inventdat = inventory.loc[inventory['USAF']==usafval]
    statyears = inventdat.loc[(inventdat['YEAR'] >= ystart) & (inventdat['YEAR'] <= yend)]
    statyears = statyears['YEAR']  

    for iyear in range(ystart,yend+1,1):
        if (iyear in statyears.astype(int).values):   
            datesofyear=[str(iyear)+'-'+i for i in datestrings]
            yearinvent = inventdat.loc[inventdat['YEAR'].astype(int) == iyear]
            yearinvent = yearinvent.loc[yearinvent['WBAN']==wbanval]
            

            fname=datpath+str(iyear)+"/"+statname[istat]+'.csv'

            try:
                data = pd.read_csv(fname)
                date_data = data[['DATE','TEMP','DEWP','STP']]
                
                # remove Feb 29th
                date_data = date_data[~date_data['DATE'].isin([str(iyear)+'-02-29'])]
                
                alldates = [str(iyear)+'-'+i for i in datestrings]
                # assign indices to dates
                alldatesinds=dict() 
                for i, j in enumerate(alldates):
                    alldatesinds.setdefault(j, []).append(i)
                    
                # find indices of all dates that are in data
                # and assign the relevant elements of t2m to the right place in the array
                res = [alldatesinds.get(i, [None]) for i in date_data['DATE']]
                t2m[istat,(iyear-ystart)*365+np.array(res).squeeze()] = date_data['TEMP']
                dewp[istat, (iyear-ystart)*365 + np.array(res).squeeze()] = date_data['DEWP']
                stp[istat, (iyear-ystart)*365 + np.array(res).squeeze()] = date_data['STP']
               
            except:
                t2m[istat,(iyear-ystart)*365:(iyear-ystart+1)*365] = nan
                dewp[istat,(iyear-ystart)*365:(iyear-ystart+1)*365] = nan
                stp[istat,(iyear-ystart)*365:(iyear-ystart+1)*365] = nan


        else:
            t2m[istat,(iyear-ystart)*365:(iyear-ystart+1)*365] = nan
            dewp[istat,(iyear-ystart)*365:(iyear-ystart+1)*365] = nan
            stp[istat,(iyear-ystart)*365:(iyear-ystart+1)*365] = nan
            

t2mxr = xr.DataArray(t2m, coords=[statname, timeout], dims=['station','time'], name='t2m')
dewpxr = xr.DataArray(dewp, coords=[statname, timeout], dims=['station','time'], name='dewp')
stpxr = xr.DataArray(stp, coords=[statname, timeout], dims=['station','time'], name='stp')
lon = xr.DataArray(lonstation, name='lon', coords=[statname], dims=['station'])
lat = xr.DataArray(latstation, name='lat', coords=[statname], dims=['station'])
stationdat = xr.merge([t2mxr, dewpxr, stpxr, lon,lat])
stationdat.time.encoding['calendar'] = "noleap"

stationdat.to_netcdf(path=fileout)






