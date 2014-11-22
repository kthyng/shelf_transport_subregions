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

if not os.path.exists('calcs'):
    os.makedirs('calcs')

whichtime = 'interannual' # 'seasonal' or 'interannual'
whichseason = 'winter' # 'winter' or 'summer' for interannual
howplot = 'log' # 'log' or 'linear'

# loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
# grid = tracpy.inout.readgrid(loc, usebasemap=True)
grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
currents_filename = list(np.sort(glob.glob('/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_????.nc')))
if 'grid' not in locals():
    grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

if whichtime == 'seasonal':

    cmap = 'YlGn'
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

        if howplot=='log':
            lev_exp = np.linspace(np.log10(0.00001),np.log10(1.0),100)
            levs = np.power(10, lev_exp)
        elif howplot=='linear':
            levs = np.linspace(0,1,100)
        mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:]/Smax, cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)

        # Outline area
        loc_shelftransport = '/home/kthyng/projects/shelf_transport/'
        STX = np.load(loc_shelftransport + 'calcs/STXpts.npz')['STX']
        xp = STX[:,0]; yp = STX[:,1]
        ax.plot(xp, yp, 'k', lw=2, alpha=0.5)

        # Horizontal colorbar below plot
        cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
        cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
        cb.set_label('Backward transport')
        levscb = np.hstack((levs[::20],levs[-1]))
        cb.set_ticks(levscb)
        cb.set_ticklabels(["%1.4f" % lev for lev in levscb])

    np.savez(fname, S=S, x=grid['xpsi'], y=grid['ypsi'])
    fig.savefig('figures/' + whichtime + '-' + howplot + '.png', bbox_inches='tight')


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

        if howplot=='log':
            lev_exp = np.linspace(np.log10(0.00001),np.log10(1.0),100)
            levs = np.power(10, lev_exp)
        elif howplot=='linear':
            levs = np.linspace(0,1,100)
        mappable = ax.contourf(grid['xpsi'], grid['ypsi'], S[i,:,:]/Smax, cmap=cmap, levels=levs, norm=colors.LogNorm())#, extend=extend)

        # Outline area
        loc_shelftransport = '/home/kthyng/projects/shelf_transport/'
        STX = np.load(loc_shelftransport + 'calcs/STXpts.npz')['STX']
        xp = STX[:,0]; yp = STX[:,1]
        ax.plot(xp, yp, 'k', lw=2, alpha=0.5)

        # Horizontal colorbar below plot
        cax = fig.add_axes([0.25, 0.075, 0.5, 0.02]) #colorbar axes
        cb = plt.colorbar(mappable, cax=cax, orientation='horizontal')
        cb.set_label('Backward transport')
        levscb = np.hstack((levs[::20],levs[-1]))
        cb.set_ticks(levscb)
        cb.set_ticklabels(["%1.4f" % lev for lev in levscb])

    np.savez(fname, S=S, x=grid['xpsi'], y=grid['ypsi'])

    fig.savefig('figures/' + whichtime + '-' + whichseason + '-' + howplot + '.png', bbox_inches='tight')



