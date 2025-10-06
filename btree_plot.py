import matplotlib.pyplot as plt
import numpy as np
from sys import exit


def plot_cells(data,ax,nodes,height,regions,col_def):
    '''Plot the binary tree grid on an existing scatter plot.'''
    col = col_def
    al = 0.8
    # Getting colours for different levels
    cols = []
    colours = ["#DC267F","#6B3BE4","#FE6100","#648FFF","#4CDC3D","#4CDC3D"]
    for x in range(height):
        [cols.append(colours[x]) for i in range(2**x)]
    # Minimum and maximum axes limits
    max_x = max(data[:,0]) + 0.07*(max(data[:,0])-min(data[:,0]))
    max_y = max(data[:,1]) + 0.07*(max(data[:,1])-min(data[:,1]))
    min_x = min(data[:,0]) - abs(0.07*(max(data[:,0])-min(data[:,0])))
    min_y = min(data[:,1]) - abs(0.07*(max(data[:,1])-min(data[:,1])))
    # Plotting the boundaries
    ax.vlines(regions[0].xmin,ymin=regions[0].ymin,ymax=max(data[:,1]),color=col,ls="--",alpha=al,lw=1.2)
    ax.hlines(regions[0].ymin,xmin=regions[0].xmin,xmax=max(data[:,0]),color=col,ls="--",alpha=al,lw=1.2)
    for r in regions:
        # col = cols[np.where(nodes==r.xmax)[0][0]] if len(np.where(nodes==r.xmax)[0])==1 else "k"
        ax.vlines(r.xmax,ymin=r.ymin,ymax=r.ymax,color=col,ls="--",alpha=al,lw=1.2)# if r.xmax<0.9*max(data[:,0]) else None
        # col = cols[np.where(nodes==r.ymax)[0][0]] if len(np.where(nodes==r.ymax)[0])==1 else "k"
        ax.hlines(r.ymax,xmin=r.xmin,xmax=r.xmax,color=col,ls="--",alpha=al,lw=1.2)# if r.ymax<0.9*max(data[:,1]) else None
        # ax.text(r.xmin+(r.xmax-r.xmin)/2,r.ymin+(r.ymax-r.ymin)/2,str(count))
    # ax.set_xlim(min_x,max_x)
    # ax.set_ylim(min_y,max_y)
    ax.invert_yaxis()


def plot_scatter_hist(x, y, obs):
    '''Plot a scatter plot with histograms on the x and y axes.'''
    def scatter_hist(x, y, ax, ax_histx, ax_histy):
        # no labels
        ax_histx.tick_params(axis="both", labelbottom=False, labelleft=False, length=0.1)
        ax_histy.tick_params(axis="both", labelbottom=False, labelleft=False, length=0.1)
        ax_histx.tick_params(axis="y", length=0)
        ax_histy.tick_params(axis="x", length=0)

        # the scatter plot:
        ax.scatter(x, y,s=25,alpha=0.25,c="grey",edgecolor=None,label="Model")

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
        lim = (int(xymax/binwidth) + 1) * binwidth

        bins = np.arange(-lim, lim + binwidth, binwidth)
        bins = 40
        ax_histx.hist(x, bins=bins,color="grey",alpha=0.5,density=True)
        ax_histy.hist(y, bins=32,orientation='horizontal',alpha=0.5,color="grey",density=True)
        # Observed data
        ax_histx.hist(obs[:,0],bins=40,color="#c9000a",density=True,alpha=0.5)
        ax_histy.hist(obs[:,1],orientation='horizontal',bins=32,color="#c9000a",density=True,alpha=0.5)

    # Start with a square Figure.
    fig = plt.figure(figsize=(8, 6))
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal axes and the main axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(2, 2,  width_ratios=(4.7, 1), height_ratios=(1, 3.7),
                        left=0.1, right=0.9, bottom=0.1, top=0.9,
                        wspace=0.05, hspace=0.062)
    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    # Delete this to get ticks back
    ax.tick_params(axis="both",length=0)

    # Limits
# Plotting 
    max_x = max(x) + 0.085*(max(x)-min(x))
    max_y = max(y) + 0.085*(max(y)-min(y))
    min_x = min(x) - abs(0.085*(max(x)-min(x)))
    min_y = min(y) - abs(0.085*(max(y)-min(y)))
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_xlabel("Separation (arcsec)")
    ax.set_ylabel(r"$\delta$mag",labelpad=-2.5)
    # Draw the scatter plot and marginals.
    scatter_hist(x, y, ax, ax_histx, ax_histy)
    return ax,ax_histx, ax_histy