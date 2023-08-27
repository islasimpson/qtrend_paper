# Routines for reading in data
import xarray as xr
import pandas as pd
import numpy as np
from pandas import Timedelta as timedelta
import sys


def fixcesmtime(dat,timebndsvar='time_bnds'):
    """ Fix the CESM timestamp using the average of time_bnds"""

    try:
        timebndavg = np.array(dat.isel(M=0)[timebndsvar],
                  dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
        dat['time'] = timebndavg
    except:
        timebndavg = np.array(dat[timebndsvar],
                  dtype='datetime64[s]').view('i8').mean(axis=1).astype('datetime64[s]')
        dat['time'] = timebndavg

    return dat


def read_sfc(filepath, datestart, dateend):
    """Read in a time slice of a surface field from datestart to dateend.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """

    #First try opening and doing the select assuming everything is working ok with the time axis
    try:
        dat = \
        xr.open_mfdataset\
        (filepath, coords="minimal", join="override", decode_times=True, use_cftime=True, compat='override').\
        sel(time=slice(datestart, dateend))
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
            print("changing longitude --> lon, latitude --> lat")
        except: pass

    except:
        #dat = xr.open_mfdataset(filepath, coords="minimal", join="override", compat='override', decode_times = False)
        dat = xr.open_mfdataset(filepath, combine='nested', concat_dim="time", coords="minimal", join="override", compat="override", decode_times = False)
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
        except: pass

        dat = xr.decode_cf(dat, use_cftime = True)
        dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    return dat
