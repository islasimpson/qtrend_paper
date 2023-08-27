import xarray as xr
import numpy as np
import sys


def compute_slope(y):
    """
    function to compute slopes along a dimension of an xarray DataArray
    combine to apply_ufunc to implement
    
    Example usage: 
    xr.apply_ufunc(linfit.compute_slope, da, vectorize=True, input_core_dims=[['time']])
    """
    x = np.arange(len(y))
    return np.polyfit(x, y, 1)[0]

def linfit_xy(x,y, sigma=None):
    """Calculate a weighted least squares linear fit between 1d arrays x and y
    Input: x (the predictor)
           y (the predictand)
           sigma (optional) the standard deviation of the y values
    Output: a and b and the residuals
    """

    if sigma is not None:
        w = 1./(sigma)
        coefs = np.polyfit(x,y,1,w=w)
    else:
        coefs = np.polyfit(x,y,1)

    a = coefs[1]
    b = coefs[0]
    return a, b

def lineardetrend(dat, dim):
    """ 
    function to detrend along a dimension of an xarray dataarray
    """

    params = dat.polyfit(dim=dim, deg=1)
    fit = xr.polyval(dat[dim], params.polyfit_coefficients)
    dat = dat-fit
    return dat


