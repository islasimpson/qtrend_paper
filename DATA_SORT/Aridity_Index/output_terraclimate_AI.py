import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan

pet_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PET_Terraclim_1980_2020_native.nc")
pet_tc = pet_tc.mean('time').load()
ppt_tc = xr.open_dataset("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/PPT_Terraclim_1980_2020_native.nc")
ppt_tc = ppt_tc.mean('time').load()
aridity = ppt_tc.ppt / pet_tc.pet
aridity = aridity.load()

aridity = aridity.rename('aridity')
aridity.to_netcdf("/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/Aridity_Index/Airidity_Index_Terraclim_native.nc")
