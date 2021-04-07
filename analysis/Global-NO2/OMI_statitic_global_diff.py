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

    period_tag1 = 'lock2'
    period_tag2 = 'lock2'
    locations = ['USA']#'Russia','Germany','USA','China','Brazil','SA','India'
    period1 = 'history'
    period2 = '2020'
    period3 = 'detrend_2020'
    tag = 'anomaly'

    path = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/periond-thre1/'
    #path_detrend = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/periond-thre1/'
    path_anomaly = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/global diff/'
    outpath = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/global diff/change/'
    trend = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/periond-thre1/detrend/'

    for country in locations:
        if tag == 'detrend':
            name1 = 'OMI-NO2-' + period1 + '-' + period_tag1 + '_' + country + '.nc4'
            name2 = 'OMI-NO2-' + period2 + '-' + period_tag1 + '_' + country + '.nc4'
            figure1 = trend + name1
            figure2 = path + name2
            NO2_1 = np.array(extract_h4_by_name(figure1, 'trend'))
            NO2_2 = np.array(extract_h4_by_name(figure2, 'time_range_Tro_NO2'))

            lon = np.array(extract_h4_by_name(figure1, 'nlon'))
            lat = np.array(extract_h4_by_name(figure1, 'nlat'))
            detrend = NO2_2 - NO2_1

            outfile = path  + 'detrend'+period2 + '_' + period_tag1 + '-' + country + '.nc4'
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
            nc_var = fid.createVariable('time_range_Tro_NO2', 'f4', ('latitude', 'longitude',))
            nc_var[:] = detrend
            nc_var.long_name = "NO2 detrend "
            nc_var.units = "10**15 molec/cm2"

            print ('finish detrend')
            fid.close()

        if tag == 'anomaly':
            name1 = 'OMI-NO2-' + period1 + '-' + period_tag1 + '_' + country + '.nc4'
            #name2 = 'OMI-NO2-mean_all'+period2+'-'+period_tag1+'_'+country+'.nc4'
            name2 = 'OMI-NO2-' + period2 + '-' + period_tag1 + '_' + country + '.nc4'
            figure1 = path + name1
            figure2 = path + name2
            NO2_1 = np.array(extract_h4_by_name(figure1, 'time_range_Tro_NO2'))
            #NO2_2 = np.array(extract_h4_by_name(figure2, 'time_range_Tro_NO2'))
            NO2_2 = np.array(extract_h4_by_name(figure2, 'time_range_Tro_NO2'))
            outfile = path_anomaly + 'anomaly_' + period2 + '_' + period1 + '_' + period_tag1 + '-' + country + '.nc4'

            lon = np.array(extract_h4_by_name(figure1, 'nlon'))
            lat = np.array(extract_h4_by_name(figure1, 'nlat'))
            difference = NO2_2 - NO2_1

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
            nc_var.long_name = "NO2 anomaly "
            nc_var.units = "10**15 molec/cm2"

            print ('finish anomaly')
            fid.close()
        if tag == 'diff':
            name1 = 'anomaly_'+period2+'_'+period1+'_'+period_tag1+'-'+country+'.nc4'
            name2 = 'anomaly_'+period2+'_'+period1+'_'+period_tag2+'-'+country+'.nc4'
            figure1 = path_anomaly + name1
            figure2 = path_anomaly + name2
            NO2_1 = np.array(extract_h4_by_name(figure1, 'NO2 anomaly'))
            NO2_2 = np.array(extract_h4_by_name(figure2, 'NO2 anomaly'))
            outfile = outpath + 'difference_of_anomaly_' +period2+'_'+ period_tag2 + '_' + period_tag1 + '-' + country + '.nc4'

            lon = np.array(extract_h4_by_name(figure1, 'nlon'))
            lat = np.array(extract_h4_by_name(figure1, 'nlat'))
            difference = NO2_2 - NO2_1

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

            nc_var = fid.createVariable('NO2 change', 'f4', ('latitude', 'longitude',))
            nc_var[:] = difference
            nc_var.long_name = "NO2 change "
            nc_var.units = "10**15 molec/cm2"

            print ('finish diff')
            fid.close()
        if tag == 'diff of change':
            name1 = 'difference_of_anomaly_'+period2+'_'+period_tag2+'_'+period_tag1+'-'+country+'.nc4'
            name2 = 'difference_of_anomaly_'+period3+'_'+period_tag2+'_'+period_tag1+'-'+country+'.nc4'
            figure1 = outpath + name1
            figure2 = outpath + name2
            NO2_1 = np.array(extract_h4_by_name(figure1, 'NO2 change'))
            NO2_2 = np.array(extract_h4_by_name(figure2, 'NO2 change'))
            outfile = outpath + 'change_' +period3+'_'+period2+'_'+ period_tag2 + '_' + period_tag1 + '-' + country + '.nc4'

            lon = np.array(extract_h4_by_name(figure1, 'nlon'))
            lat = np.array(extract_h4_by_name(figure1, 'nlat'))
            difference = NO2_2 - NO2_1

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

            nc_var = fid.createVariable('NO2 change', 'f4', ('latitude', 'longitude',))
            nc_var[:] = difference
            nc_var.long_name = "NO2 change "
            nc_var.units = "10**15 molec/cm2"

            print ('finish diff of change')
            fid.close()


