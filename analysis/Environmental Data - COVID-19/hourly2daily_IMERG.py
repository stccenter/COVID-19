'''
Author: your name
Date: 2020-09-24 15:17:10
LastEditTime: 2020-10-08 15:20:33
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: \Environmental Data - COVID-19\hourly2daily_IMERG.py
'''
import glob, os,sys
import numpy as np
import h5py
from scipy.io import netcdf
from netCDF4 import Dataset



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


    infolder='D:/data process/COVID-19/analysis/Environmental Data - COVID-19/in_forder/'
    output_dir='D:/data process/COVID-19/analysis/Environmental Data - COVID-19/out_forder/'
    #all_nc_files, n_all=get_files(infolder,'*.nc4')
    #all_nc_files=np.sort(all_nc_files)
    # Put all the desired date into date_range
    date_range = ['20200601']

    for idate in date_range:
        all_nc_files, n_all = get_files(infolder, '3B-HHR-E.MS.MRG.3IMERG.'+idate+'*')
        allday_prec = []
        print(n_all)
        print(idate)
        print(infolder)
        #print (all_nc_files)
        for i in range(n_all):
            the_filename = all_nc_files[i]

            thePrec = np.array(extract_h4_by_name(the_filename, 'Grid/precipitationCal')).squeeze()
            thePrec = np.transpose(thePrec)
            # print(thePrec)
            # exit()
            #print (thePrec.shape)
            #exit()
            lats = extract_h4_by_name(the_filename, 'Grid/lat')
            lons = extract_h4_by_name(the_filename, 'Grid/lon')
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




