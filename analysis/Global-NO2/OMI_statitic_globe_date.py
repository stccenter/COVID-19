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
from lunardate import LunarDate

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

    periods = ['2020']
    flags = ['after2','lock2']#'before',,'after'
    locations = ['USA']#'Russia','Germany','China','Brazil','SA','India'
    lock_year = 2020
    output_dir = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/period/'
    for location in locations:
        if location == 'China':
            lock_month = 1
            lock_day = 25
            reopen_month = 4
            reopen_day = 9
            end_month = 6
            end_day = 30
            lock_doy = doy(lock_year, lock_month, lock_day)
            lunar_lock_date = LunarDate.fromSolarDate(int(lock_year), int(lock_month), int(lock_day))
            lstart = LunarDate(lunar_lock_date.year, 1, 1)
            ldoy_lock = (lunar_lock_date - lstart).days + 1

            reopen_doy = doy(lock_year, reopen_month, reopen_day)
            lunar_reopen_date = LunarDate.fromSolarDate(int(lock_year), int(reopen_month), int(reopen_day))
            ldoy_reopen = (lunar_reopen_date - lstart).days + 1
            print (ldoy_reopen)

            end_doy = doy(lock_year, end_month, end_day)
            lunar_end_date = LunarDate.fromSolarDate(int(lock_year), int(end_month), int(end_day))
            ldoy_end = (lunar_end_date - lstart).days + 1

            start = 330
            lockdown = 355
            reopen = ldoy_reopen
            end = ldoy_end
        if location == 'USA':
            start_month = 1
            start_day = 1
            lock_month = 3
            lock_day = 21
            reopen_month = 6
            reopen_day = 14
            lock2_month = 11
            lock2_day = 17
            end_month = 12
            end_day = 30
            start = doy(lock_year, start_month, start_day)
            lockdown = doy(lock_year, lock_month, lock_day)
            reopen = doy(lock_year, reopen_month, lock2_day)
            lock2 = doy(lock_year, lock2_month, reopen_day)
            end = doy(lock_year, end_month, end_day)
        if location == 'Brazil':
            start_month = 1
            start_day = 1
            lock_month = 3
            lock_day = 21
            reopen_month = 6
            reopen_day = 1
            end_month = 6
            end_day = 30
            start = doy(lock_year, start_month, start_day)
            lockdown = doy(lock_year, lock_month, lock_day)
            reopen = doy(lock_year, reopen_month, reopen_day)
            end = doy(lock_year, end_month, end_day)
        if location == 'SA':
            start_month = 1
            start_day = 1
            lock_month = 3
            lock_day = 26
            reopen_month = 6
            reopen_day = 7
            end_month = 6
            end_day = 30
            start = doy(lock_year, start_month, start_day)
            lockdown = doy(lock_year, lock_month, lock_day)
            reopen = doy(lock_year, reopen_month, reopen_day)
            end = doy(lock_year, end_month, end_day)
        if location == 'India':
            start_month = 1
            start_day = 1
            lock_month = 3
            lock_day = 22
            reopen_month = 5
            reopen_day = 4
            end_month = 6
            end_day = 30
            start = doy(lock_year, start_month, start_day)
            lockdown = doy(lock_year, lock_month, lock_day)
            reopen = doy(lock_year, reopen_month, reopen_day)
            end = doy(lock_year, end_month, end_day)
        if location == 'Germany':
            start_month = 1
            start_day = 1
            lock_month = 3
            lock_day = 22
            reopen_month = 5
            reopen_day = 19
            end_month = 6
            end_day = 30
            start = doy(lock_year, start_month, start_day)
            lockdown = doy(lock_year, lock_month, lock_day)
            reopen = doy(lock_year, reopen_month, reopen_day)
            end = doy(lock_year, end_month, end_day)
        if location == 'Russia':
            start_month = 1
            start_day = 1
            lock_month = 3
            lock_day = 19
            reopen_month = 6
            reopen_day = 1
            end_month = 6
            end_day = 30
            start = doy(lock_year, start_month, start_day)
            lockdown = doy(lock_year, lock_month, lock_day)
            reopen = doy(lock_year, reopen_month, reopen_day)
            end = doy(lock_year, end_month, end_day)

        for period in periods:
            for flag in flags:
                if location == 'China':
                    infolder='/storage/qliu6/trmm/1998_2017_merge/results/OMI/'+period+'_daily/'+period+'-lunardate-thre1/'
                else:
                    infolder = '/storage/qliu6/trmm/1998_2017_merge/results/OMI/' + period + '_daily/' + period + '-solardate-thre1/'
                all_nc_files, n_all=get_files(infolder,'*.nc4')
                all_nc_files=np.sort(all_nc_files)
                all_NO2 = []
                for i in range(n_all):
                    the_filename = all_nc_files[i]
                    idoy = int(the_filename.split('.nc4')[0][-3:])
                    if flag == 'before':
                        if idoy not in range(start,lockdown):
                            continue
                    if flag == 'during':
                        if location == 'China':
                            if idoy not in range(1, reopen):
                                continue
                        else:
                            if idoy not in range(lockdown, reopen):
                                continue
                    if flag == 'after':
                        if idoy not in range(reopen,lock2):
                            continue
                    if flag == 'lock2':
                        if idoy not in range(lock2, end + 1):
                            continue
                    theQV2M = np.array(extract_h4_by_name(the_filename, 'Tro_NO2'))
                    idx_invalid = np.where(theQV2M<0)
                    theQV2M[idx_invalid] = np.nan
                    theQV2M = np.array(theQV2M)
                    theQV2M = theQV2M
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
                mean_NO2 = np.nanmean(all_NO2,axis=0)
                N = all_NO2.shape[0]
                mean_NO2 = np.nanmean(all_NO2, axis=0)
                stand_error = np.nanstd(all_NO2, axis=0) / math.sqrt(N)

                outfile = output_dir+ 'OMI-NO2-'+period+'-'+flag+'_'+location+'.nc4'
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

                nc_var = fid.createVariable('time_range_Tro_NO2', 'f4', ('latitude','longitude',))
                nc_var[:] = mean_NO2
                nc_var.long_name = "average tropospheric NO2"
                nc_var.units = "10**15 molec/cm2"

                nc_var = fid.createVariable('time_range_Tro_NO2_standard_error', 'f4', ('latitude', 'longitude',))
                nc_var[:] = stand_error
                nc_var.long_name = "standard error average tropospheric NO2"
                nc_var.units = "10**15 molec/cm2"

                fid.close()
                print ('finish ...')

    '''
    period = '2020'
    flag = 'before'
    location = 'China'
    infolder='/storage/qliu6/trmm/1998_2017_merge/results/OMI/'+period+'_daily/'+period+'-lunardate/'
    infolder_so2 = '/storage/qliu6/OMI/history/acdisc.gesdisc.eosdis.nasa.gov/data/Aura_OMI_Level3/OMSO2e.003/' + period + '/'
    output_dir='/storage/qliu6/trmm/1998_2017_merge/results/OMI/period/'
    all_nc_files, n_all=get_files(infolder,'*.nc4')
    all_nc_files=np.sort(all_nc_files)
    months = ['01','02','03','04','05','06']


    lock_year = 2020
    lock_month = 1
    lock_day = 25
    lock_doy = doy(lock_year,lock_month,lock_day)
    reopen_month = 4
    reopen_day = 9
    reopen_doy = doy(lock_year,reopen_month,reopen_day)
    end_month = 6
    end_day = 30
    end_doy = doy(lock_year, end_month, end_day)

    all_NO2 = []
    for i in range(n_all):
        the_filename = all_nc_files[i]

        #date = the_filename.split('_'+period+'m')[1][:4]
        #imonth = date[:2]
        #iday = date[2:4]
        #idoy = doy(int(period),int(imonth),int(iday))
        idoy = int(the_filename.split('.nc')[0][-3:])
        if flag == 'before':
            if idoy not in range(0,lock_doy):
                continue
        if flag == 'during':
            if idoy not in range(lock_doy, reopen_doy):
                continue
        if flag == 'after':
            if idoy not in range(reopen_doy, end_doy):
                continue
        print (the_filename)
        #theQV2M = np.array(extract_h4_by_name(the_filename, '/HDFEOS/GRIDS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop'))
        theQV2M = np.array(extract_h4_by_name(the_filename, 'Tro_NO2'))
        idx_invalid = np.where(theQV2M<0)
        theQV2M[idx_invalid] = np.nan
        theQV2M = np.array(theQV2M)
        theQV2M = theQV2M/math.pow(10,15)
        #lats = extract_h4_by_name(the_filename, '/HDFEOS/GRIDS/lat')
        #lons = extract_h4_by_name(the_filename, '/HDFEOS/GRIDS/lon')

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
    mean_NO2 = np.nanmean(all_NO2,axis=0)
    N = all_NO2.shape[0]
    mean_NO2 = np.nanmean(all_NO2, axis=0)
    stand_error = np.nanstd(all_NO2, axis=0) / math.sqrt(N)


    outfile = output_dir+ 'OMI-NO2-'+period+'-'+flag+'_'+location+'.nc4'
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

    nc_var = fid.createVariable('time_range_Tro_NO2', 'f4', ('latitude','longitude',))
    nc_var[:] = mean_NO2
    nc_var.long_name = "average tropospheric NO2"
    nc_var.units = "10**15 molec/cm2"

    nc_var = fid.createVariable('time_range_Tro_NO2_standard_error', 'f4', ('latitude', 'longitude',))
    nc_var[:] = stand_error
    nc_var.long_name = "standard error average tropospheric NO2"
    nc_var.units = "10**15 molec/cm2"


    fid.close()


    print ('finish ...')
    '''



