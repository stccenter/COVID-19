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


    infolder='/datavolume/covid19/precipitation/hourly/'
    output_dir='/datavolume/covid19/precipitation/daily/'
    #all_nc_files, n_all=get_files(infolder,'*.nc4')
    #all_nc_files=np.sort(all_nc_files)
    # Put all the desired date into date_range
    date_range = ['20200102','20200103','20200104','20200105','20200106',
          '20200107','20200108','20200109','20200110','20200111','20200112',
          '20200113','20200114','20200115','20200116','20200117','20200118',
          '20200119','20200120','20200121','20200122','20200123','20200124',
          '20200125','20200126','20200127','20200128','20200129','20200130','20200131',
        '20200201','20200202','20200203','20200204','20200205','20200206',
                  '20200207','20200208','20200209','20200210','20200211','20200212',
                  '20200213','20200214','20200215','20200216','20200217','20200218',
                  '20200219','20200220','20200221','20200222','20200223','20200224',
                  '20200225','20200226','20200227''20200228','20200229',
        '20200301','20200302','20200303','20200304','20200305','20200306',
                  '20200307','20200308','20200309','20200310','20200311','20200312',
                  '20200313','20200314','20200315','20200316','20200317','20200318',
                  '20200319','20200320','20200321','20200322','20200323','20200324',
                  '20200325','20200326','20200327','20200328','20200329','20200330','20200331',
                  '20200401','20200402','20200403','20200404','20200405','20200406',
          '20200407','20200408'
         ]

    for idate in date_range:
        all_nc_files, n_all = get_files(infolder, '3B-HHR-E.MS.MRG.3IMERG.'+idate+'*')
        allday_prec = []
        #print (all_nc_files)
        for i in range(n_all):
            the_filename = all_nc_files[i]

            thePrec = np.array(extract_h4_by_name(the_filename, 'precipitationCal')).squeeze()
            thePrec = np.transpose(thePrec)
            #print (thePrec.shape)
            #exit()
            lats = extract_h4_by_name(the_filename, 'lat')
            lons = extract_h4_by_name(the_filename, 'lon')
            isif = thePrec.shape
            allday_prec.append(thePrec)
        allday_prec = np.array(allday_prec)

        daily_prec = np.nanmean(allday_prec, axis=0)

        outfile = output_dir+ 'daily_precipitation_'+idate+'.nc4'
        # create nc file
        fid = netcdf.netcdf_file(outfile, 'w')
        # create dimension variable, so we can use it in the netcdf
        fid.createDimension('longitude', isif[1])
        fid.createDimension('latitude', isif[0])

        nc_var = fid.createVariable('nlat', 'f4', ('latitude',))
        nc_var[:] = lats
        nc_var.long_name = 'latitude'
        nc_var.standard_name = 'latitude'
        nc_var.units = 'degrees_north'

        nc_var = fid.createVariable('nlon', 'f4', ('longitude',))
        nc_var[:] = lons
        nc_var.long_name = 'longitude'
        nc_var.standard_name = 'longitude'
        nc_var.units = 'degrees_east'

        nc_var = fid.createVariable('daily_precipitation', 'f4', ('latitude','longitude',))
        nc_var[:] = daily_prec
        nc_var.long_name = "daily_precipitation"
        nc_var.units = "mm/hr"


        fid.close()


        print ('finish ...')
        #exit()




