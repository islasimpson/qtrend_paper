import importlib
import xarray as xr
import numpy as np
import xesmf as xe
import sys
import warnings

# DATA for interpolation grid
cesmdat = xr.open_dataset(
    '../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/"
wgtfile=pathout+'wgtfile.nc'

path="/project/mojave/observations/TerraClimate/pet/"
pet = xr.open_mfdataset(path+"*.nc")
pet = pet.sel(time=slice("1980-01-01","2020-12-31"))

pet = pet.mean('time')

#regridder = xe.Regridder(pet, grid_out, 'bilinear', periodic=True, reuse_weights=False,
#filename=wgtfile)

#pet_rg = regridder(pet)
pet.to_netcdf(pathout+'PET_Terraclim_1980_2020_native.nc')

