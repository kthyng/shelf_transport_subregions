''' Script to run drifters at pre-selected subregions of the domain to
get stream function plots out. '''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
from datetime import datetime, timedelta
import glob
from tracpy.tracpy_class import Tracpy
from matplotlib.mlab import find, Path


# loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'
grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
currents_filename = list(np.sort(glob.glob('/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_????.nc')))
grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

def init(name):
    '''
    Initialization for the simulation.
    '''

    # loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'

    time_units = 'seconds since 1970-01-01'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 25 # in-between tracks: 12 # old tracks: 25 

    # Number of steps to divide model output for outputting drifter location
    N = 1 # not interested in drifter locations

    # Number of days
    ndays = 30

    # This is a backward-moving simulation
    ff = -1 

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0. # old tracks: 5.
    av = 0. # m^2/s

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines.
    dostream = 1

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid_filename=grid_filename, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, dostream=dostream, savell=False, doperiodic=0, 
                N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, 
                time_units=time_units, usebasemap=True, grid=grid, vert_filename=vert_filename)

    # tp._readgrid()

    # Initial lon/lat locations for drifters
    startindsfile = 'CH-starting-inds.npz'
    if not os.path.exists(startindsfile):
        # get starting drifter locations for this region 
        loc_shelftransport = '/home/kthyng/projects/shelf_transport/'
        CH = np.load(loc_shelftransport + 'calcs/CHpts.npz')['CH']
        xp = CH[:,0]; yp = CH[:,1]
        gpath = Path(np.vstack((xp, yp)).T)
        drifters = netCDF.Dataset(loc_shelftransport + 'tracks/2004-01-01T00gc.nc')
        xg = drifters.variables['xg'][:]; yg = drifters.variables['yg'][:]
        # Change to projected drifter locations now
        nanind = np.isnan(xg) + (xg==-1) # indices where nans are location in xg, yg; for reinstitution of nans
        # pdb.set_trace()
        xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2xy') 
        xp[nanind] = np.nan; yp[nanind] = np.nan
        del(xg,yg) # don't need grid info anymore
        # save indices of drifters that start in the coastal areas
        inds = gpath.contains_points(np.vstack((xp[:,0].flat, yp[:,0].flat)).T).reshape(xp[:,0].shape)
        np.savez(startindsfile, inds=inds)
    else:
        inds = np.load(startindsfile)['inds']

    startptsfile = 'CH-starting-ll.npz'
    if not os.path.exists(startptsfile):
        loc_shelftransport = '/home/kthyng/projects/shelf_transport/'
        drifters = netCDF.Dataset(loc_shelftransport + 'tracks/2004-01-01T00gc.nc')
        xg = drifters.variables['xg'][:]; yg = drifters.variables['yg'][:]
        xg = xg[inds,0]; yg = yg[inds,0]
        lon0, lat0, _ = tracpy.tools.interpolate2d(xg, yg, tp.grid, 'm_ij2ll')
        np.savez(startptsfile, lon0=lon0, lat0=lat0)
    else:
        d = np.load(startptsfile)
        lon0 = d['lon0']; lat0 = d['lat0']
        d.close()

    # equal weightings for drifters for transport.
    T0 = np.ones(lon0.size, order='F')

    U = np.ma.zeros(tp.grid['xu'].shape, order='F')
    V = np.ma.zeros(tp.grid['xv'].shape, order='F')

    # pdb.set_trace()
       
    return tp, lon0, lat0, T0, U, V


def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')
        
    #loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    #grid = tracpy.inout.readgrid(loc)

    overallstartdate = datetime(2009, 7, 1, 0, 1)
    overallstopdate = datetime(2009, 9, 1, 0, 1)

    date = overallstartdate

    # Start from the beginning and add days on for loop
    # keep running until we hit the next month
    while date < overallstopdate:

        name = date.isoformat()[0:13]

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc') and \
            not os.path.exists('tracks/' + name + 'gc.nc'):

            # Read in simulation initialization
            tp, lon0, lat0, T0, U, V = init(name)

            # Run tracpy
            # Save directly to grid coordinates
            lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0, T0=T0, U=U, V=V)


        # Increment by 24 hours for next loop, to move through more quickly
        date = date + timedelta(hours=24)

    # # Do transport plot
    # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=startdate.isoformat()[0:7] + '*', 
    #     extraname=startdate.isoformat()[0:7], Title='Transport on Shelf', dmax=1.0)


if __name__ == "__main__":
    run()    

