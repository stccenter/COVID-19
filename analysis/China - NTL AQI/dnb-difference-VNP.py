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
#import geopandas as gpd
#import pysal as ps
#from pysal.contrib.viz import mapping as maps
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
    #plg_dir = '/storage/qliu6/VIIRS/covid19/Hubei_province/'
    #plg_file = plg_dir + 'CHN_admin2_hubei.shp'
    #lsoas = gpd.read_file(plg_file)
    #print (lsoas.head())
    #exit()

    #input_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/nc/'
    #output_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/diff/'
    input_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/nc/'
    output_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/diff/'
    #pre_nc_filename = 'China202001_mean.nc4'
    #post_nc_filename = 'China202002_mean.nc4'
    region = 'China'
    #region = 'Wuhan'
    pre_nc_filename = 'China201812_mean.nc4'
    post_nc_filename = 'China201902_mean_SFD.nc4'
    #date1 = pre_nc_filename.split('.')[0][-8:]
    #date2 = post_nc_filename.split('.')[0][-8:]
    date1 = pre_nc_filename.split('_m')[0][-6:]
    date2 = post_nc_filename.split('_m')[0][-6:]

    pre_file = input_dir + pre_nc_filename
    post_file = input_dir + post_nc_filename

    pre_rad =np.array(extract_h4_by_name(pre_file, 'monthly_mean_radiance'))
    idx_bg1 = np.where(pre_rad<5)
    pre_rad[idx_bg1] = 0
    print (np.nanmean(pre_rad))
    post_rad =np.array(extract_h4_by_name(post_file, 'monthly_mean_radiance'))
    idx_bg2 = np.where(post_rad<5)
    post_rad[idx_bg2] = 0
    print (np.nanmean(post_rad))
    #exit()
    lats =np.array(extract_h4_by_name(pre_file, 'nlat'))
    lons = np.array(extract_h4_by_name(pre_file, 'nlon'))
    
    rad_diff = post_rad - pre_rad
    sif = rad_diff.shape
    n_pixels = sif[0]

    '''''
    outfile_nc = output_dir + region + '_'+date2+' and '+date1+'.nc4'
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

    nc_var = fid.createVariable('radiance difference', 'f4', ('n_pixels',))
    nc_var[:] = rad_diff
    nc_var.units = 'nW/(cm2 sr)'

    fid.close()
    '''''
    #exit()

    outfile = output_dir  + region + '_'+date2+' and '+date1+'.png'
    y_min = lats.min()
    y_max = lats.max()
    x_min = lons.min()
    x_max = lons.max()
    y_min = 30.18
    y_max = 31.12
    x_min = 113.78
    x_max = 114.9
    ratio = (x_max - x_min) / (y_max - y_min)
    mindata = rad_diff[~np.isnan(rad_diff)].min()
    maxdata = rad_diff[~np.isnan(rad_diff)].max()


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
    plt.title(region + '_' + 'Radiance difference between '+ date2+' and '+date1)
    # convert lat and lon to map projection coordinates
    lons, lats = map(lons, lats)

    # create render
    cmap = mpl.cm.seismic  # rainbow
    maxdata = 10  # NY
    mindata = -10
    ticker_width = maxdata - mindata
    nticks = 10
    ticker_interval = ticker_width / nticks
    nticks += 1
    normticks = np.arange(mindata, maxdata, ticker_interval)
    norm = mpl.colors.Normalize(vmin=mindata, vmax=maxdata, clip=True)

    cs = map.scatter(lons, lats, s=5, marker='s', c=rad_diff, cmap=cmap, norm=norm, edgecolors='none')
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
    map.drawmeridians(np.arange(0., 420., x_ivl), linewidth=0.0, color='k', labels=[False, False, False, True])
    # set x title and y title
    ax.set_xlabel("", fontsize=5)
    ax.set_ylabel("", fontsize=5)
    # draw boundary lines
    # map.drawstates(color='black', linewidth=1, ax=ax)
    # map.drawcountries(color='black', linewidth=0.5, ax=ax)
    # map.drawcoastlines(linewidth=1, color='black', ax=ax)
    # map.readshapefile('/storage/qliu6/VIIRS/covid19/Hubei_province/', 'comarques')
    cbar.set_label('Radiance nW/(cm2 sr)' + '\n' + '\n' + 'Data Provided by NOAA CLASS, NASA NCCS and GMAO' +
                   '\n' + 'Visual Analytics Conducted by Qian Liu, NSF Spatiotemporal Innovation Center')

    plt.savefig(outfile, bbox_inches='tight', dpi=300)
    print ('finish')

