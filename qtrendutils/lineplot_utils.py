import numpy as np
import matplotlib.pyplot as plt

def plotlinetime_j2d_monthly(fig, data, x1, x2, y1, y2, titlestr, yrange=None,
    yticks=None, yticklabels=None, ytitle=None, linecolor=None, points=True, label=None, linestyle='solid'):
    """ plot a line plot using monthly data from Jan to Dec
        Input: fig = your figure 
           data = a 365 element array containing data to be plotted
           x1 = location of left edge of plot
           x2 = location of right edge of plot
           y1 = location of bottom edge of plot
           y2 = location of top edge of plot
           titlestr = plot title
           yrange = optional range for y axis
           yticks = optional ticks for y axis
           yticklabels = optional tick labels for y axis
           ytitle= optional title for y axis
           linecolor = optional color of line
    """

    ax = fig.add_axes([x1,y1,x2-x1,y2-y1])

    monticks = np.arange(0,12,1)
    monticks2 = np.arange(0,12,1)+0.5
    if (yrange):
        ax.set_ylim(yrange)

    if (yticks):
        ax.set_yticks(yticks)

    if (yticklabels):
        ax.set_yticklabels(yticklabels, fontsize=13)

    if (ytitle):
        ax.set_ylabel(ytitle, fontsize=14)

    ax.set_xlim([0,12])
    ax.tick_params(which='minor', length=0)
    ax.set_xticks(monticks)
    ax.set_xticklabels([])
    ax.set_xticks(monticks2, minor=True)
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], minor=True, fontsize=13)
    ax.set_title(titlestr, fontsize=16)

    xvals = np.arange(0,12,1)+0.5
    xpad = np.zeros([len(data)+2]).astype('float')
    xpad[0] = xvals[len(xvals)-1]-12 ; xpad[len(xpad)-1] = xvals[0]+12 ; xpad[1:len(xpad)-1] = xvals
    datpad = np.zeros([len(data)+2]).astype('float')
    datpad[0] = data[len(data)-1].values 
    datpad[len(xpad)-1] = data[0].values 
    datpad[1:len(xpad)-1] = data.values

    if linecolor is not None:
            ax.plot(xpad,datpad,color=linecolor, linewidth=2, label=label, linestyle=linestyle)
            if (points == True):
                ax.plot(xpad,datpad,"o",markerfacecolor=linecolor,
                markeredgecolor="black", markersize=10)
    else:
        ax.plot(xpad,datpad, linewidth=2, label=label, linestyle=linestyle)
        if (points == True):
            ax.plot(xpad,datpad,"o",markeredgecolor="black",
            markersize=10, markeredgewidth=2)

    return ax

def oplotlinetime_j2d_monthly(ax, data, linecolor=None, linewidth=1, points=True, label=None):
    """ overplot a line on a January - December monthly line plot"""

    xvals = np.arange(0,12,1)+0.5
    xpad = np.zeros([len(data)+2]).astype('float')
    xpad[0] = xvals[len(xvals)-1]-12 ; xpad[len(xpad)-1] = xvals[0]+12 ; xpad[1:len(xpad)-1] = xvals

    datpad = np.zeros([len(data)+2]).astype('float')
    datpad[0] = data[len(data)-1].values 
    datpad[len(xpad)-1] = data[0].values
    datpad[1:len(xpad)-1] = data.values

    if (linecolor[0]):
        ax.plot(xpad,datpad,color=linecolor, linewidth=linewidth, label=label)
        if (points == True):
            ax.plot(xpad,datpad,"o",markerfacecolor=linecolor,
             markeredgecolor="black", markersize=10)
    else:
        ax.plot(xpad,datpad, linewidth=linewidth, label=label)
        if (points == True):
            ax.plot(xpad,datpad,"o",markeredgecolor="black",
            markersize=10, markeredgewidth=2)

    return ax

