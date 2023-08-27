import numpy as np
import xarray as xr
import sys

def runningmean(dat, nysm, timeaxis='time', dropna=False):
    """dat = your data with a time axis with name equal to whatever you set "timeaxis" to
       nysm = the number of time values in your running mean
       dropna = False if you don't want to drop the NaN's at the edges
    """

    window_kwargs = {timeaxis:nysm}
    if (dropna):
        datm = dat.rolling(center=True, min_periods=nysm, **window_kwargs).mean(timeaxis).dropna(timeaxis)
    else:
        datm = dat.rolling(center=True, min_periods=nysm, **window_kwargs).mean(timeaxis)
    return datm

