'''
Plot some drifter tracks and location histograms, as determined by input indices.
'''

import numpy as np
import matplotlib.pyplot as plt
import glob
import netCDF4 as netCDF
import tracpy
import tracpy.plotting
import matplotlib as mpl
import pdb
import op
from matplotlib import ticker, colors, cm
from matplotlib.mlab import find
import os


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


# plottracks = False # plot drifter trajectories
# plothist = True # plot drifter locations at a time as a histogram

whichtime = 'interannual' # 'seasonal' or 'interannual'
whichseason = 'winter' # 'winter' or 'summer' for interannual

# # which drifters to plot for summer and winter, respectively
# inds = np.load('calcs/summer-drifter-indices.npz')['inds']
# indw = np.load('calcs/winter-drifter-indices.npz')['inds']

loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc, usebasemap=True)

# # decimate the drifter indices
# ddi = 1000

if whichtime == 'seasonal':

    cmap = 'YlGn'

    shelf_depth = 100

    fname = 'calcs/' + whichtime + '-S.npz'

    S = np.zeros((2,grid['xpsi'].shape[0],grid['xpsi'].shape[1])) # initialize

    fig, axarr = plt.subplots(1,2)
    fig.set_size_inches(13.675, 6.6125)
    fig.subplots_adjust(left=0.04, bottom=0.15, right=1.0, top=0.96, wspace=0.07, hspace=0.04)

    for i, ax in enumerate(axarr):

        # Titles for subplots
        if i==0:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
            ax.set_title('Winter')
            Files = glob.glob('tracks/20??-0[1,2]-*gc.nc')
        elif i==1:
            tracpy.plotting.background(grid=grid, ax=ax, parslabels=[0,0,0,0], mers=np.arange(-100, -80, 2))
            ax.set_title('Summer')
            Files = glob.glob('tracks/20??-0[7,8]-*gc.nc')

        if not os.path.exists(fname):
            for File in Files:
                # print File
                d = netCDF.Dataset(File)
                # pdb.set_trace()
                U = d.variables['U'][:]; V = d.variables['V'][:]
                d.close()
                Stemp = np.sqrt(op.resize(U, 1)**2 + op.resize(V, 0)**2)
                S[i,:,:] = S[i,:,:] + Stemp

            # locator = ticker.MaxNLocator(11)
            # locator.create_dummy_axis()
            # locator.set_bounds(0, 1) 
            # levels = locator()
            # extend = 'max'
            # H = H/Hmax
            Smax = 1.

        else:
            d = np.load(fname); S = d['S']; d.close()
            Smax = S.max()

        lev_exp = np.linspace(np.log10(1000), 5.25,50)
        # lev_exp = np.linspace(np.log10(1000), np.ceil(np.log10(S.max())),50)
        # lev_exp = np.arange(1e-13, np.ceil(np.log10(S.max())+1))
        levs = np.power(10, lev_exp)

        mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:], cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)
        ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)
        # pdb.set_trace()

        # outline the area where drifters started
        d = np.load('../../shelf_transport/calcs/winter-contour-pts.npz')
        ax.plot(d['x'], d['y'], 'k', lw=3)
        d.close()

    np.savez(fname, S=S, x=grid['xpsi'], y=grid['ypsi'])

    fig.savefig('figures/' + whichtime + '.png', bbox_inches='tight')


elif whichtime == 'interannual':

    cmap = 'YlGn'

    shelf_depth = 100

    fname = 'calcs/' + whichtime + '-' + whichseason + '-S.npz'

    S = np.zeros((7,grid['xpsi'].shape[0],grid['xpsi'].shape[1])) # initialize

    fig, axarr = plt.subplots(2,4)
    fig.set_size_inches(13.4, 6.6125)
    fig.subplots_adjust(left=0.03, bottom=0.15, right=1.0, top=0.96, wspace=0.03, hspace=0.11)

    for i, ax in enumerate(axarr.flatten()):

        yr = 2004+i
        if whichseason == 'winter':
            Files = glob.glob('tracks/' + str(yr) + '-0[1,2]-*gc.nc')
        elif whichseason == 'summer':
            Files = glob.glob('tracks/' + str(yr) + '-0[7,8]-*gc.nc')

        # Titles for subplots
        if i==4:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2))
        elif i==7:
            ax.set_frame_on(False)
            ax.set_axis_off()
            continue
        else:
            tracpy.plotting.background(grid=grid, ax=ax, mers=np.arange(-100, -80, 2), 
                merslabels=[0, 0, 0, 0], parslabels=[0, 0, 0, 0])

        ax.set_title(str(yr))

        if not os.path.exists(fname):
            for File in Files:
                # print File
                d = netCDF.Dataset(File)
                # pdb.set_trace()
                U = d.variables['U'][:]; V = d.variables['V'][:]
                d.close()
                Stemp = np.sqrt(op.resize(U, 1)**2 + op.resize(V, 0)**2)
                S[i,:,:] = S[i,:,:] + Stemp

            # locator = ticker.MaxNLocator(11)
            # locator.create_dummy_axis()
            # locator.set_bounds(0, 1) 
            # levels = locator()
            # extend = 'max'
            # H = H/Hmax
            Smax = 1.

        else:
            d = np.load(fname); S = d['S']; d.close()
            Smax = S.max()

        lev_exp = np.linspace(np.log10(1000), 4.5,50)
        # lev_exp = np.linspace(np.log10(1000), np.ceil(np.log10(S.max())),50)
        # lev_exp = np.arange(1e-13, np.ceil(np.log10(S.max())+1))
        levs = np.power(10, lev_exp)
        # pdb.set_trace()
        mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:], cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)
        # mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:]/Smax, cmap=cmap, levels=np.linspace(0,0.25), extend='max')
        ax.contour(grid['xr'], grid['yr'], grid['h'], [shelf_depth], colors='0.1', linewidth=3)

        # outline the area where drifters started
        d = np.load('../../shelf_transport/calcs/winter-contour-pts.npz')
        ax.plot(d['x'], d['y'], 'k', lw=3)
        d.close()

    np.savez(fname, S=S, x=grid['xpsi'], y=grid['ypsi'])

    fig.savefig('figures/' + whichtime + whichseason + '.png', bbox_inches='tight')



