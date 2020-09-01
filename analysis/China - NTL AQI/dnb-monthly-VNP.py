import glob, os,sys
# import netCDF4
import math
import scipy.interpolate
import scipy.ndimage
import numpy as np
import matplotlib as mpl
import scipy.stats as stats
import matplotlib.mlab as mlab
import h5py
from scipy.io import netcdf
from scipy.stats import lognorm
from scipy.stats import gamma
from scipy.stats import chisquare
#from compiler.ast import flatten
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
from scipy.stats import norm
from netCDF4 import Dataset
import pandas as pd
from sklearn.decomposition import PCA
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap,addcyclic, shiftgrid,cm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import matplotlib
import datetime
import h5py
matplotlib.use('Agg')



def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
     
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.
 
    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates
 
    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin
 
    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)
 
    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print ("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = np.asarray( newdims, dtype=int )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa
 
    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]
 
        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )
 
        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )
 
            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )
 
        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )
 
        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = np.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = np.mgrid[nslices]
 
        newcoords_dims = range(np.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords
 
        newcoords_tr += ofs
 
        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas
 
        newcoords_tr -= ofs
 
        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print ("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None

def get_files(dir,ext):
    allfiles=[]
    os.chdir(dir)
    for file in glob.glob(ext):
        allfiles.append(file)
    return allfiles, len(allfiles)

    # read netcdf 3 file by dataset name
def extract_nc3_by_name(filename, dsname):
    nc_data = netcdf.netcdf_file(filename, "r")
    ds = np.array(nc_data.variables[dsname][:])
    nc_data.close()
    return ds
def extract_h5_by_name(filename,dsname):
    h5_data = h5py.File(filename)
    ds = np.array(h5_data[dsname][:])
    h5_data.close()
    return ds

def extract_h4_by_name(filename,dsname):
    h4_data = Dataset(filename)
    ds = np.array(h4_data[dsname][:])
    h4_data.close()
    return ds
# -----------------------------------------------------------------------------
# -||||||||||||||||||||||||Main function|||||||||||||||||||||||||||||||||||||||
# -----------------------------------------------------------------------------
if __name__ == '__main__':

    input_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/nc/'
    output_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/png/'
    alldate = [
        '20200101','20200102','20200103','20200104','20200105','20200106',
'20200107','20200108','20200109','20200110','20200111','20200112',
'20200113','20200114','20200115','20200116','20200117','20200118',
'20200119','20200120','20200121','20200122','20200123','20200131',
               ]
    month_name = 'January'
    region = 'China'
    #region = 'Wuhan'
    #region = 'Hubei'
    #region = 'nan'
    all_radiance = []

    for date in alldate:
        year = int(date[:4])
        month = int(date[4:6])
        the_file = input_dir+'night_rad_'+ region+date+'.nc'
        all_nc_files, n_all = get_files(input_dir, the_file)
        if n_all ==0:
            continue
        theRadiance = extract_h4_by_name(the_file, 'Radiance')
        lats = extract_h4_by_name(the_file, 'nlat')
        lons = extract_h4_by_name(the_file, 'nlon')
        cloud_mask = extract_h4_by_name(the_file, 'Cloud mask')
        isif = theRadiance.shape
        idx_cloud = np.where(cloud_mask > 64)
        theRadiance[idx_cloud]=np.nan
        #idx_nan = np.where(theRadiance==np.nan)
        print (date)
        all_radiance.append(theRadiance)

    all_radiance = np.array(all_radiance)
    all_radiance_mean = np.nanmean(all_radiance, axis=0)
    sif = all_radiance_mean.shape
    print (sif)
    n_pixels = sif[0]

    #idx_background = np.where(all_radiance_mean < 5)
    #all_radiance_mean[idx_background] = 0
    all_radiance_mean = np.array(all_radiance_mean)


    outfile_nc = input_dir + region + str(year)+ str('{0:02}'.format(month)) + '_mean_SFD.nc4'
    # create nc file
    fid = netcdf.netcdf_file(outfile_nc, 'w')
    # create dimension variable, so we can use it in the netcdf
    fid.createDimension('n_pixels', n_pixels)

    nc_var = fid.createVariable('nlat', 'f4', ('n_pixels',))
    nc_var[:] = lats
    nc_var.long_name = 'latitude'
    nc_var.standard_name = 'latitude'
    nc_var.units = 'degrees_north'

    nc_var = fid.createVariable('nlon', 'f4', ('n_pixels',))
    nc_var[:] = lons
    nc_var.long_name = 'longitude'
    nc_var.standard_name = 'longitude'
    nc_var.units = 'degrees_east'

    nc_var = fid.createVariable('monthly_mean_radiance', 'f4', ('n_pixels',))
    nc_var[:] = all_radiance_mean
    nc_var.units = 'nW/(cm2 sr)'

    fid.close()
    exit()

    outfile = output_dir  + region + '_'+ str(year) + str('{0:02}'.format(month)) + 'hubei.png'
    y_min = lats.min()
    y_max = lats.max()
    x_min = lons.min()
    x_max = lons.max()
    # lon_west = 113.78
    # lon_east = 114.9
    # lat_south = 30.18
    # lat_north = 31.12
    #y_min = 30.18
    #y_max = 31.12
    #x_min = 113.78
    #x_max = 114.9
    y_min = 29
    y_max = 33.38
    x_min = 108
    x_max = 116.38

    ratio = (x_max - x_min) / (y_max - y_min)
    idx_nan = np.where(all_radiance_mean==np.nan)
    #all_radiance_mean[idx_nan] = 0
    all_radiance_mean = np.array(all_radiance_mean)
    mindata = all_radiance_mean[~np.isnan(all_radiance_mean)].min()
    maxdata = all_radiance_mean[~np.isnan(all_radiance_mean)].max()


    # print(x_min)
    # print(x_max)
    # print(y_min)
    # print(y_max)


    fig = plt.figure(figsize=(16, 12))  # Create a new figure window
    rect = [0.125, 0.25, 0.5, 0.5 / ratio]  # [left, bottom, width, height] (ratio 0~1)0.256
    ax = plt.axes(rect)

    # create a basemap
    map = Basemap(projection='cyl', llcrnrlat=y_min, urcrnrlat=y_max, \
                  llcrnrlon=x_min, urcrnrlon=x_max, ax=ax)  # lon_0=0.0,
    plt.title('Hubei' + '_' + 'Monthly Mean Night Light Radiance of '+ str(year) + ' ' + month_name)
    # convert lat and lon to map projection coordinates
    lons, lats = map(lons, lats)

    # create render
    cmap = mpl.cm.jet # seismic
    #cmap = mpl.cm.gist_gray
    maxdata = 60
    mindata = 0
    ticker_width = maxdata - mindata
    nticks = 10
    ticker_interval = ticker_width / nticks
    nticks += 1
    normticks = np.arange(mindata, maxdata, ticker_interval)
    norm = mpl.colors.Normalize(vmin=mindata, vmax=maxdata, clip=True)

    cs = map.scatter(lons, lats, s=5, marker='s', c=all_radiance_mean, cmap=cmap, norm=norm, edgecolors='none')
    # create color bar
    cbaxes = fig.add_axes([rect[0], rect[1] - 0.05, rect[2], 0.02])
    # divider = make_axes_locatable(ax)
    # cbaxes = divider.append_axes("bottom", size="5%", pad=0.05)
    cbar = fig.colorbar(cs, cmap=cmap, ax=ax, cax=cbaxes, orientation='horizontal',
                        ticks=normticks, fraction=0.046, pad=0.05, extend='both', extendfrac='auto')  #
    # cbar.set_label('Some Units')
    tick_locator = mpl.ticker.MaxNLocator(nbins=nticks)
    cbar.locator = tick_locator
    cbar.ax.xaxis.set_ticks_position('bottom')

    # draw grid
    y_ivl = (x_max - x_min) / 5
    map.drawparallels(np.arange(-90., 90., y_ivl), linewidth=0.0, color='k', labels=[True, False, False, False])
    x_ivl = (y_max - y_min) / 5
    map.drawmeridians(np.arange(-180., 180., x_ivl), linewidth=0.0, color='k', labels=[False, False, False, True])
    # set x title and y title
    ax.set_xlabel("", fontsize=3)
    ax.set_ylabel("", fontsize=3)
    # draw boundary lines
    # map.drawstates(color='black', linewidth=1, ax=ax)
    # map.drawcountries(color='black', linewidth=0.5, ax=ax)
    # map.drawcoastlines(linewidth=1, color='black', ax=ax)
    # map.readshapefile('/storage/qliu6/VIIRS/covid19/Hubei_province/', 'comarques')
    cbar.set_label('Radiance nW/(cm2 sr)' + '\n' + '\n' + 'Data Provided by NOAA CLASS, NASA NCCS and GMAO' +
                   '\n' + 'Visual Analytics Conducted by Qian Liu, NSF Spatiotemporal Innovation Center')

    plt.savefig(outfile, bbox_inches='tight', dpi=300)
    print ('finish')

