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

def get_subfiles(dir, ext):
    allfiles = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if (file.endswith(ext)):
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

def is_leap_year(year):
    """ if year is a leap year return True
        else return False """
    if year % 100 == 0:
        return year % 400 == 0
    return year % 4 == 0

def doy(Y,M,D):
    """ given year, month, day return day of year
        Astronomical Algorithms, Jean Meeus, 2d ed, 1998, chap 7 """
    if is_leap_year(Y):
        K = 1
    else:
        K = 2
    N = int((275 * M) / 9.0) - K * int((M + 9) / 12.0) + D - 30
    return N

def ymd(Y,N):
    """ given year = Y and day of year = N, return year, month, day
        Astronomical Algorithms, Jean Meeus, 2d ed, 1998, chap 7 """
    if is_leap_year(Y):
        K = 1
    else:
        K = 2
    M = int((9 * (K + N)) / 275.0 + 0.98)
    if N < 32:
        M = 1
    D = N - int((275 * M) / 9.0) + K * int((M + 9) / 12.0) + 30
    return Y, M, D

# -----------------------------------------------------------------------------
# -||||||||||||||||||||||||Main function|||||||||||||||||||||||||||||||||||||||
# -----------------------------------------------------------------------------
if __name__ == '__main__':

    period = '2020'
    years = [2020]#2010,2011,2012,2013,2014,2105,2016,2017,2018,
    infolder = '/storage/qliu6/OMI/history/acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/OMNO2d.003/'
    output_dir = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/'+period+'_daily/'+period+'-solardate/'
    annual_file = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/OMI-NO2-global-annualmean.nc4'
    annual_mean = np.array(extract_h4_by_name(annual_file, 'monthly_Tro_NO2'))
    all_nc_files, n_all = get_subfiles(infolder, '.he5')
    all_nc_files = np.sort(all_nc_files)
    doys = np.linspace(183,365,183)

    for j in doys:
        all_NO2 = []
        for i in range(n_all):
            the_filename = all_nc_files[i]
            year = the_filename.split('OMNO2d_')[1][:4]
            month = the_filename.split('OMNO2d_')[1][5:7]
            day = the_filename.split('OMNO2d_')[1][7:9]
            if int(year) not in years:
                continue
            idoy = doy(int(year), int(month), int(day))
            if idoy != j:
                continue
            theQV2M = np.array(extract_h4_by_name(infolder+year+'/'+the_filename, '/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop'))
            idx_nan = np.where(annual_mean < 1)
            theQV2M[idx_nan] = np.nan
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
            all_NO2.append(theQV2M)
        all_NO2 = np.array(all_NO2)
        N = all_NO2.shape[0]
        X = np.linspace(1, N, N)
        n_lat = all_NO2.shape[1]
        n_lon = all_NO2.shape[2]
        trend = np.zeros((n_lat, n_lon))
        for m in range(n_lat):
            for n in range(n_lon):
                the_TS = np.array(all_NO2[:, m, n])
                not_nan_ind = ~np.isnan(the_TS)
                if np.array(the_TS[not_nan_ind]).size == 0:
                    trend[m, n] = np.nan
                    continue
                else:
                    TS_valid = np.array(the_TS[not_nan_ind])
                    # the_TS_detrend = signal.detrend(TS_valid)
                    slope, b, r_val, p_val, std_err = stats.linregress(X[not_nan_ind], TS_valid)
                    trend[m, n] = slope

        mean_NO2 = np.nanmean(all_NO2, axis=0)
        stand_error = np.nanstd(all_NO2,axis=0)/math.sqrt(N)

        outfile = output_dir+ 'OMI-NO2-'+'{0:03}'.format(int(j))+'.nc4'
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

        nc_var = fid.createVariable('Tro_NO2', 'f4', ('latitude','longitude',))
        nc_var[:] = mean_NO2
        nc_var.long_name = "average tropospheric NO2"
        nc_var.units = "10**15 molec/cm2"

        nc_var = fid.createVariable('Tro_NO2_standard_error', 'f4', ('latitude', 'longitude',))
        nc_var[:] = stand_error
        nc_var.long_name = "standard error average tropospheric NO2"
        nc_var.units = "10**15 molec/cm2"
        fid.close()
        print ('finish '+ str(j))
        
        #exit()




