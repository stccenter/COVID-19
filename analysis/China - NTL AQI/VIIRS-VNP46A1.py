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
import shapefile
import matplotlib
import datetime
from pyhdf.SD import SD, SDC
from pyhdf.HDF import *
from pyhdf.VS import *
import h5py
matplotlib.use('Agg')

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

def read_vpn46a1_boundary(file):
    with h5py.File(file, mode='r') as f:
        group = f['/HDFEOS/GRIDS/VNP_Grid_DNB']
        min_lon = group.attrs['WestBoundingCoord'][0]
        max_lon = group.attrs['EastBoundingCoord'][0]
        min_lat = group.attrs['SouthBoundingCoord'][0]
        max_lat = group.attrs['NorthBoundingCoord'][0]
        return np.array([min_lon,max_lon,min_lat,max_lat])
# -----------------------------------------------------------------------------
# -||||||||||||||||||||||||Main function|||||||||||||||||||||||||||||||||||||||
# -----------------------------------------------------------------------------
if __name__ == '__main__':

    #Data_infolder = '/storage/qliu6/VIIRS/covid19/DNB500/NY/storage/qliu6/VIIRS/covid19/DNB500/NY/'
    #Data_infolder = '/home/qliu6/storage/qliu6/VIIRS/covid19/DNB500/WH/'
    Data_infolder = '/storage/qliu6/VIIRS/covid19/China-gridded/'
    output_dir = '/storage/qliu6/VIIRS/covid19/results/VNP46A1/'
    alldate = [ '20190301','20190302','20190303','20190304','20190305','20190306',
'20190307','20190308','20190309','20190310','20190311','20190312',
'20190313','20190314','20190315','20190316','20190317','20190318',
'20190319','20190320','20190321','20190322','20190323','20190324',
'20190325','20190326','20190327','20190328','20190329','20190330','20190331',
               ]

    #region = 'NY'
    #lon_west = -74.5
    #lon_east = -71.8
    #lat_south = 40.4
    #lat_north = 41.2

    #region = 'Wuhan'
    #lon_west = 113.68
    #lon_east = 115.1
    #lat_south = 29.96
    #lat_north = 31.37

    #lon_west = 113.78
    #lon_east = 114.9
    #lat_south = 30.18
    #lat_north = 31.12

    #region = 'nan'
    #lon_west = 120
    #lon_east = 126
    #lat_south = 34
    #lat_north = 40

    region = 'China'
    lon_west = 73
    lon_east = 135
    lat_south = 18
    lat_north = 53.52
    ratio = (lon_east - lon_west) / (lat_north - lat_south)

    for date in alldate:
        year = int(date[:4])
        month = int(date[4:6])
        day = int(date[6:8])
        doy = (datetime.datetime(year, month, day)-datetime.datetime(year,1,1)).days+1
        doy = '{0:03}'.format(doy)
        all_nc_files, n_all = get_files(Data_infolder, 'VNP46A1.A' + str(year)+str(doy)+'*.h5')
        all_nc_files = np.sort (all_nc_files)
        if n_all !=35:
            continue
        j = 0
        for i in range(n_all):

            the_filename = all_nc_files[i]
            #print (the_filename)
            the_file = Data_infolder + the_filename
            theRadiance = 0.1*extract_h5_by_name(the_file, '/HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/DNB_At_Sensor_Radiance_500m')
            isif = theRadiance.shape
            n_lon = isif[0]
            n_lat = isif[1]
            coors = read_vpn46a1_boundary(the_file)
            max_lon = coors[1]
            min_lon = coors[0]
            max_lat = coors[3]
            min_lat = coors[2]

            x = np.linspace(min_lon, max_lon,n_lon)
            y = np.linspace(max_lat,min_lat,n_lat)
            lon_grid,lat_grid = np.meshgrid(x,y)

            theRadiance = theRadiance.reshape(-1)
            idx_nan=np.where(theRadiance<0)
            theRadiance[idx_nan]=np.nan
            lats = lat_grid.reshape(-1)
            lons = lon_grid.reshape(-1)
            zenith = 0.01*extract_h5_by_name(the_file, 'HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/Solar_Zenith')
            zenith = zenith.reshape(-1)
            cloud_mask = extract_h5_by_name(the_file, 'HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/QF_Cloud_Mask')
            cloud_mask = cloud_mask.reshape(-1)
            cloud_mask = cloud_mask & 0b00011000000

            region_lon = np.logical_and(lons > lon_west, lons < lon_east)
            region_lat = np.logical_and(lats > lat_south, lats < lat_north)
            region_loc = np.logical_and(region_lon, region_lat)
            idx_cloud = np.where(cloud_mask>64)
            #theRadiance[idx_cloud]=np.nan
            theRadiance = np.array(theRadiance)
            #idx_clear = np.logical_and(cloud_mask<=64, region_loc)
            #idx_region, = np.where(idx_clear)
            idx_region, = np.where(region_loc)

            if len(idx_region)==0:
                continue

            theRadiance = theRadiance[idx_region]
            lats = lats[idx_region]
            lons = lons[idx_region]
            cloud_mask = cloud_mask[idx_region]

            theRadiance = np.array(theRadiance)
            lats = np.array(lats)
            lons = np.array(lons)
            cloud_mask = np.array(cloud_mask)
            isif = theRadiance.shape
            #print (isif)

            if j == 0:
                global_radiance = theRadiance
                global_lat = lats
                global_lon = lons
                global_cloud = cloud_mask
                j = 1
            else:
                global_radiance = np.concatenate((global_radiance, theRadiance), axis=0)
                global_lat = np.concatenate((global_lat, lats), axis=0)
                global_lon = np.concatenate((global_lon, lons), axis=0)
                global_cloud = np.concatenate((global_cloud, cloud_mask), axis=0)

        global_radiance = np.array(global_radiance)
        global_lat = np.array(global_lat)
        global_lon = np.array(global_lon)
        global_cloud = np.array(global_cloud)
        sif = global_radiance.shape
        # print (sif)
        if sif[0] == 0:
            continue
        print (sif)
        #save figure
        y_min=global_lat.min()
        y_max=global_lat.max()
        x_min=global_lon.min()
        x_max=global_lon.max()
        mindata=global_radiance[~np.isnan(global_radiance)].min()
        maxdata=global_radiance[~np.isnan(global_radiance)].max()

            # save as nc file
        outfile_nc = output_dir + 'nc/night_rad_' + region + date + '.nc'
        fid = netcdf.netcdf_file(outfile_nc, 'w')
        # create dimension variable, so we can use it in the netcdf
        fid.createDimension('npixels', sif[0])

        # latitude
        nc_var = fid.createVariable('nlat', 'f4', ('npixels',))
        nc_var[:] = global_lat
        nc_var.long_name = 'latitude'
        nc_var.standard_name = 'latitude'
        nc_var.units = 'degrees_north'
        # longitude
        nc_var = fid.createVariable('nlon', 'f4', ('npixels',))
        nc_var[:] = global_lon
        nc_var.long_name = 'longitude'
        nc_var.standard_name = 'longitude'
        nc_var.units = 'degrees_east'

        # radiance
        nc_var = fid.createVariable('Radiance', 'f4', ('npixels',))
        nc_var[:] = global_radiance
        nc_var.units = 'nW/(cm2 sr)'

        nc_var = fid.createVariable('Cloud mask', 'f4', ('npixels',))
        nc_var[:] = global_cloud
        nc_var.units = ''
        # end output
        fid.close()
        print ('finish '+ date)





