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
import math


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

#To calculate difference between two time periods
    period1 = 'peri-peri-hist'
    period2 = 'post-post-hist'


    path = '/home/qliu6/trmm/1998_2017_merge/results/OMI/'
    name1 = 'difference_'+period1+'.nc4'
    name2 = 'difference_'+period2+'.nc4'
    figure1 = path + name1
    figure2 = path + name2
    NO2_1 = np.array(extract_h4_by_name(figure1, period1+'_Tro_NO2'))
    NO2_2 = np.array(extract_h4_by_name(figure2, period2+'_Tro_NO2'))
    lon = np.array(extract_h4_by_name(figure1, 'nlon'))
    lat = np.array(extract_h4_by_name(figure1, 'nlat'))
    difference = NO2_2 - NO2_1
    #ratio = difference/NO2_1

    outfile = path + 'anomaly_' + period2+'-'+period1 + '.nc4'
    # create nc file
    fid = netcdf.netcdf_file(outfile, 'w')
    # create dimension variable, so we can use it in the netcdf
    fid.createDimension('longitude', lon.shape[0])
    fid.createDimension('latitude', lat.shape[0])

    nc_var = fid.createVariable('nlat', 'f4', ('latitude',))
    nc_var[:] = lat
    nc_var.long_name = 'latitude'
    nc_var.standard_name = 'latitude'
    nc_var.units = 'degrees_north'

    nc_var = fid.createVariable('nlon', 'f4', ('longitude',))
    nc_var[:] = lon
    nc_var.long_name = 'longitude'
    nc_var.standard_name = 'longitude'
    nc_var.units = 'degrees_east'

    nc_var = fid.createVariable('NO2 anomaly', 'f4', ('latitude', 'longitude',))
    nc_var[:] = difference
    nc_var.long_name = "NO2 difference between "+period2+" and "+period1
    nc_var.units = "10**15 molec/cm2"

    fid.close()

    print ('finish ...')
    exit()


#To calculate periodic mean

    period = 'post'
    infolder='/storage/qliu6/OMI/'+period+'/'
    output_dir='/home/qliu6/trmm/1998_2017_merge/results/OMI/'
    all_nc_files, n_all=get_files(infolder,'*.he5')
    all_nc_files=np.sort(all_nc_files)

    all_NO2 = []
    for i in range(n_all):
        the_filename = all_nc_files[i]
        theQV2M = np.array(extract_h4_by_name(the_filename, '/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop'))
        idx_invalid = np.where(theQV2M<0)
        theQV2M[idx_invalid] = np.nan
        theQV2M = np.array(theQV2M)
        theQV2M = theQV2M/math.pow(10,15)
        min_lon = -180.0
        max_lon = 180.0
        max_lat = -90.0
        min_lat = 90.0
        n_lon = 1440
        n_lat = 720
        x = np.linspace(min_lon, max_lon, n_lon)
        y = np.linspace(max_lat, min_lat, n_lat)


        #print (theQV2M.shape,x.shape)
        all_NO2.append(theQV2M)

    all_NO2 = np.array(all_NO2)
    mean_NO2 = np.nanmean(all_NO2,axis=0)
    #print (mean_NO2.shape)

    outfile = output_dir+ 'OMI-NO2-'+period+'.nc4'
    # create nc file
    fid = netcdf.netcdf_file(outfile, 'w')
    # create dimension variable, so we can use it in the netcdf
    fid.createDimension('longitude', x.shape[0])
    fid.createDimension('latitude', y.shape[0])

    nc_var = fid.createVariable('nlat', 'f4', ('latitude',))
    nc_var[:] = y
    nc_var.long_name = 'latitude'
    nc_var.standard_name = 'latitude'
    nc_var.units = 'degrees_north'

    nc_var = fid.createVariable('nlon', 'f4', ('longitude',))
    nc_var[:] = x
    nc_var.long_name = 'longitude'
    nc_var.standard_name = 'longitude'
    nc_var.units = 'degrees_east'

    nc_var = fid.createVariable(period+'_Tro_NO2', 'f4', ('latitude','longitude',))
    nc_var[:] = mean_NO2
    nc_var.long_name = "average tropospheric NO2"
    nc_var.units = "10**15 molec/cm2"


    fid.close()


    print ('finish ...')




