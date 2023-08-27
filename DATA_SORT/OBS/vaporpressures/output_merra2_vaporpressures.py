import xarray as xr
import numpy as np
import xesmf as xe

from qtrendutils import readdata_utils as read
from qtrendutils import humiditycalcs as qcalcs

import warnings
warnings.filterwarnings('ignore')

pathout="/project/cas/islas/python_savs/qtrend_paper/DATA_SORT/OBS/vaporpressures/"
tfile="/project/mojave/observations/MERRA2/mon/T2M_MERRA2_asm_mon_198001_202012.nc"
qfile="/project/mojave/observations/MERRA2/mon/QV2M_MERRA2_asm_mon_198001_202012.nc"
psfile="/project/mojave/observations/MERRA2/mon/PS_MERRA2_asm_mon_198001_202012.nc"

tdat_native = read.read_sfc(tfile, "1980-01","2020-12")
qdat_native = read.read_sfc(qfile, "1980-01","2020-12")
psdat_native = read.read_sfc(psfile,"1980-01","2020-12")

cesmdat = xr.open_dataset('../../cesm_grid.nc')
grid_out = xr.Dataset({'lat': (['lat'], cesmdat.lat.values)}, {'lon': (['lon'], cesmdat.lon.values)})

wgtfile = pathout+"wgtfile.nc"

regridder = xe.Regridder(tdat_native, grid_out, 'bilinear', periodic=True, reuse_weights=False,
                 filename=wgtfile)
tdat = regridder(tdat_native.T2M)
qdat = regridder(qdat_native.QV2M)
psdat = regridder(psdat_native.PS)

svp = qcalcs.calcsvp(tdat)
svp = svp.rename('svp')

vp = qcalcs.calcvpfromhuss(qdat, psdat)
vp = vp.rename('vp')

vp.to_netcdf(path=pathout+'vaporpressures_MERRA2.nc')
svp.to_netcdf(path=pathout+'vaporpressures_MERRA2.nc', mode='a')
