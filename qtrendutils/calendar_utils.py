## routines for calculating seasonal climatology and seasonal timeseries
import xarray as xr
import numpy as np
from datetime import timedelta, datetime
import pandas as pd
from math import nan
import sys

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

def leap_year(year, calendar='standard'):
    """Determine if year is a leap year
    Args: 
        year (numeric)
    """
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

def get_days_per_mon(time, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    
    Args: time (CFTimeIndex): ie. ds.time.to_index()
          calendar (str): default 'standard'
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if ( (leap_year(year, calendar=calendar)) and (month == 2)):
            month_length[i] += 1
    return month_length


def YYYYMM2date(time, caltype='standard'):
    """ Convert a date of the form YYYYMM to a datetime64 object """
    date = pd.to_datetime(time, format='%Y%m')
    return date

def YYYYMMDD2date(date, caltype='standard'):
    time = pd.to_datetime(date, format='%Y%m%d')
    return time


def calcannualmean(ds, skipna=False):
    """ Calculate the annual mean weighting the months of the year appropriately if
        given the calendar type
    """

    def dothecalc(var, skipna=False):
        month_length = var.time.dt.days_in_month
        wghts = month_length.groupby('time.year') / month_length.groupby('time.year').sum()
        if (skipna):
            datsum = (var*wghts).groupby('time.year').sum(dim='time', skipna=True)
            cond = var.isnull()
            ones = xr.where(cond, 0, 1)
            onesum = (ones*wghts).groupby('time.year').sum(dim='time')
        else:
            datsum = (var*wghts).groupby('time.year').sum(dim='time', skipna=False)
            cond = var.isnull()
            ones = xr.where(cond, 1, 1)
            onesum = (ones*wghts).groupby('time.year').sum(dim='time')

        var_am = datsum / onesum 
        return var_am

    #--Note if the NaN's are different in each variable you'll be averaging over
    #-- different times for each variable
    dset = False
    try:
        varnames = list(ds.keys())
        dset = True
    except:
        pass

    if (dset):
        for i, ivar in enumerate(varnames):
            var = ds[ivar]
            var_am = dothecalc(var, skipna=skipna)
            var_am = var_am.rename(ivar)
            if (i == 0):
                ds_am = var_am
            else:
                ds_am = xr.merge([ds_am, var_am])
    else:
        ds_am = dothecalc(ds, skipna=skipna)


    return ds_am 


def group_monthly2yearly(dat):
    years = dat['time.year']
    ybeg = np.array(years[0])
    yend = np.array(years[len(years)-1])

    datyear=[]
    for iyear in np.arange(ybeg,yend+1,1):
        yearlydat = dat.sel(time=slice(str(iyear)+'-01-01', str(iyear)+'-12-31'))
        yearlydat['time'] = np.arange(0,12,1)
        datyear.append(yearlydat)

    datyear = xr.concat(datyear, dim='year')
    datyear['year'] = np.arange(ybeg,yend+1,1) 
    return datyear


