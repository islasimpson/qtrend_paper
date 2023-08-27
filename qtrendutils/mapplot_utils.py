import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from qtrendutils import colormap_utils as mycolors

import cartopy as cart
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm
from matplotlib import colors as c

def contourmap_northamerica_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, markersize=20, signifdat=None):
    """ plot a map plot of scatter points for the northern america 
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)

    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)
    ax.set_extent([-170, -50, 10, 80], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", s=markersize, vmin=cmin, vmax=cmax, cmap = mymap)

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if ( signifdat is not None ):
      ax.scatter(lon, lat, c=signifdat, marker="o", s=markersize, vmin=0, vmax=2,
                 cmap='gray_r')
      #ax.scatter(lon, lat, c=signifdat, marker="o", color='gray')

    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax


def contourmap_southwest_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, markersize=20, signifdat=None):
    """ plot a map plot of scatter points for the northern america 
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)

    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE)
    ax.set_extent([-130, -100, 30, 43], crs = ccrs.PlateCarree())

#    if (labels):
#        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
#        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
#        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
#        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
#        ax.xformatter = LongitudeFormatter()
#        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    ax.scatter(lon, lat, c=dat, marker="o", s=markersize, vmin=cmin, vmax=cmax, cmap = mymap)

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if ( signifdat is not None ):
      ax.scatter(lon, lat, c=signifdat, marker="o", s=markersize, vmin=0, vmax=2,
                 cmap='gray_r')
      #ax.scatter(lon, lat, c=signifdat, marker="o", color='gray')

    #ax.scatter(lon, lat, c=dat, marker="o", vmin=-170, vmax=170, cmap="RdYlBu_r")
    return ax




def contourmap_northamerica_fill_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", maskocean=False, contourlines=False, contourlinescale=1):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.PlateCarree())
    ax.set_aspect('auto')
    ax.add_feature(cfeature.COASTLINE, zorder=100)
    ax.set_extent([-170, -50, 10, 80], crs = ccrs.PlateCarree())

    if (labels):
        ax.set_xticks([-150, -100, -50], crs = ccrs.PlateCarree())
        ax.set_xticklabels(['150W','100W','50W'], fontsize=12)
        ax.set_yticks([20,40,60,80], crs = ccrs.PlateCarree())
        ax.set_yticklabels(['20N','40N','60N','80N'], fontsize=12)
        ax.xformatter = LongitudeFormatter()
        ax.yformatter = LatitudeFormatter()

    ax.set_title(titlestr, fontsize=16)

    dat, lon = add_cyclic_point(dat, coord=lon)
    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap)

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if (contourlines):
        clevlines = clevs*contourlinescale
        clevlines = clevlines[np.abs(clevlines) > ci/2. ]
        ax.contour(lon,lat,dat, levels=clevlines, colors='black', transform=ccrs.PlateCarree())


    return ax


def contourmap_bothcontinents_robinson_pcolormesh_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, cmap="blue2red", fontsize=15, signifdat=None, onecolor=False, color=None, oplot=False, ax=None):
    """ plot pcolormesh on a robinson map"""

    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    norm = BoundaryNorm(clevs, ncolors=mymap.N, clip=True)

    if (oplot == False):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
        ax.set_aspect('auto')
        ax.set_title(titlestr, fontsize=fontsize)
        ax.set_extent([-180,180,-57,90], crs = ccrs.PlateCarree())

    
    if (onecolor):
        cmap = c.ListedColormap([color])
        ax.pcolormesh(lon, lat, dat, cmap=cmap, norm=norm, edgecolor='none',
           transform=ccrs.PlateCarree())
    else:
        ax.pcolormesh( lon, lat, dat, cmap=mymap, norm=norm, edgecolor='none', 
         transform=ccrs.PlateCarree())

    if ( signifdat is not None ):
        mygraymap = c.ListedColormap(['lightgray','lightgray'])
        graynorm = BoundaryNorm([0,1], ncolors=mygraymap.N, clip=True)
        ax.pcolormesh(lon,lat,signifdat,cmap=mygraymap, norm=graynorm, 
              transform=ccrs.PlateCarree())

    if (oplot == False):
        ax.add_feature(cfeature.COASTLINE, zorder=100)

    ax.axis('off')


    return ax


def contourmap_bothcontinents_robinson_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red", onecolor=False, color=None, fontsize=15, maskocean=False, oplot=False, ax=None,
 markersize=10):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs)

    if (oplot == False):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
        ax.set_aspect('auto')
        ax.add_feature(cfeature.COASTLINE)

        ax.set_title(titlestr, fontsize=fontsize)

    if (onecolor):
        ax.scatter(lon, lat, c=color, s=markersize, marker="o",
               transform=ccrs.PlateCarree(), zorder=100) 
    else:
        ax.scatter(lon, lat, c=dat, s=markersize, marker="o", vmin=cmin, vmax=cmax, cmap = mymap, 
               transform=ccrs.PlateCarree(), zorder=100)

    if (~oplot):

        ax.set_global()

        if (maskocean):
            ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    return ax


def contourmap_continentsonly_robinson_noborder_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red",nowhite=False, fontsize=15, maskocean=False, signifdat=None,
oplot=False, onecolor=False, color='black',stipplesignif=False):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)


    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs, nowhite=nowhite)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs, nowhite=nowhite)
 

    if (oplot == False):
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
        ax.set_aspect('auto')
        ax.add_feature(cfeature.COASTLINE, zorder=100.)

        ax.set_title(titlestr, fontsize=fontsize)

    dat, lon = add_cyclic_point(dat, coord=lon)
    if (onecolor):
        ax.contourf(lon, lat, dat, level=clevs, colors=color, edgecolor='none', 
             transform=ccrs.PlateCarree())
    else:
        ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", 
              transform=ccrs.PlateCarree())

    ax.set_global()

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    if (oplot == False):
        ax.axis('off')
        ax.set_extent([-180,180,-57,90], crs = ccrs.PlateCarree())

    if ( signifdat is not None ):
        lonsignif = signifdat.lon
        signifdat, lonsignif = add_cyclic_point( signifdat, coord=lonsignif)
        if (stipplesignif):
            density=3
            ax.contourf(lonsignif, lat, signifdat, levels=[0,0.5,1], colors='none',
                 hatches=[density*'.',density*'.', density*','],
               transform = ccrs.PlateCarree())
        else:
            ax.contourf(lonsignif, lat, signifdat, levels=[0,0.5,1], colors='lightgray',
            transform = ccrs.PlateCarree())


    return ax


def contourmap_continentsonly_robinson_noborder_scatter_pos(fig, dat, lon, lat, ci, cmin, cmax, titlestr,
 x1, x2, y1, y2, labels=True, cmap="blue2red",nowhite=False, fontsize=15, maskocean=False, signifdat=None,
 onecolor=False, color=None, markersize=10, oplot=False, ax=None):
    """ plot a contour map of 2D data dat with coordinates lon and lat
        Input:
              fig = the figure identifier
              dat = the data to be plotted
              lon = the longitude coordinate
              lat = the latitude coordinate
              ci = the contour interval
              cmin = the minimum of the contour range
              cmax = the maximum of the contour range
              titlestr = the title of the map
              x1 = position of the left edge
              x2 = position of the right edge
              y1 = position of the bottom edge
              y2 = position of the top edge
              labels = True/False (ticks and  labels are plotted if true) 
              cmap = color map (only set up for blue2red at the moment)
    """

#    if (signifdat is not None):
#        markersize=markersize/5

    # set up contour levels and color map
    nlevs = (cmax-cmin)/ci + 1
    clevs = np.arange(cmin, cmax+ci, ci)

    if (cmap == "blue2red"):
        mymap = mycolors.blue2red_cmap(nlevs, nowhite=nowhite)

    if (cmap == "precip"):
        mymap = mycolors.precip_cmap(nlevs, nowhite=nowhite)

    if (oplot == False): 
        ax = fig.add_axes([x1, y1, x2-x1, y2-y1], projection=ccrs.Robinson(central_longitude=0))
        ax.set_aspect('auto')
        ax.add_feature(cfeature.COASTLINE)

        ax.set_title(titlestr, fontsize=fontsize)

    if (onecolor):
        ax.scatter(lon, lat, c=color, s=markersize, marker="o",
               transform=ccrs.PlateCarree(), zorder=0)
    else:
        ax.scatter(lon, lat, c=dat, s=markersize, marker="o", vmin=cmin, vmax=cmax, cmap = mymap,
               transform=ccrs.PlateCarree(), zorder=0)

#    ax.scatter(lon, lat, c=dat, s=markersize, marker="o", vmin=cmin, vmax=cmax, cmap=mymap,
#                transform=ccrs.PlateCarree(), zorder=100)

#    dat, lon = add_cyclic_point(dat, coord=lon)
#    ax.contourf(lon, lat, dat, levels=clevs, cmap = mymap, extend="both", transform=ccrs.PlateCarree())

    ax.set_global()

    if (maskocean):
        ax.add_feature(cart.feature.OCEAN, zorder=1, edgecolor='k')

    ax.axis('off')
    ax.set_extent([-180,180,-57,90], crs = ccrs.PlateCarree())

    if ( signifdat is not None ):
        lonsignif = lon[~np.isnan(signifdat)]
        latsignif = lat[~np.isnan(signifdat)]
        datsignif = dat[~np.isnan(signifdat)]
#        print(lonsignif)

        ax.scatter(lonsignif, latsignif, c=datsignif, s=markersize, marker="o", 
               vmin=cmin, vmax=cmax, cmap = mymap,
               transform=ccrs.PlateCarree(), zorder=0)
  

#        ax.scatter(lonsignif, latsignif, facecolors='gray',s=markersize, edgecolors='gray', 
#            transform = ccrs.PlateCarree(), alpha=0.2)


    return ax



