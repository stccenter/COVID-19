
import glob, os, sys
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal, gdal_array, osr, ogr
import rasterio
from rasterio.mask import mask
from geopandas import *;
import geopandas as geopd


def get_files(dir, ext):
    allfiles = []
    os.chdir(dir)
    for file in glob.glob(ext):
        allfiles.append(file)
    return allfiles, len(allfiles)


def extract_h4_by_name(filename, dsname):
    h4_data = Dataset(filename)
    #    print(h4_data)
    ds = np.array(h4_data[dsname][:])
    h4_data.close()
    return ds


"""
search all the vectors, return the statistical 2D array
"""


def readShapefileFrom2DNCFile(shpdir, datas, xReSolution, yReSolution, ncData, keyFieldName):
    shpdata = geopd.GeoDataFrame.from_file(shpdir)
    # get the location of factors
    xyValues = [[keyFieldName, 'Max', 'Mean', 'Min', 'std']]

    for i in range(0, len(shpdata)):
        # get features of vector data
        geo = shpdata.geometry[i]
        feature = [geo.__geo_interface__]
        itemindex = np.argwhere(shpdata.columns == keyFieldName)[0][0]
        name = shpdata.iloc[i][itemindex]
        #        print(len(shpdata.columns),str(shpdata.columns[itemindex]),name,itemindex)
        # cut the .nc4 file through features, rasterio.mask.mask
        out_image, out_transform = mask(datas, feature, all_touched=True, crop=True, nodata=datas.nodata)
        # get overlap area of features
        #        out_image, out_transform = mask(datas, geometry, all_touched=True, crop=True, nodata=datas.nodata)
        values = out_image.tolist()
        values = values[0]

        # delete errors
        mindata = ncData[~np.isnan(ncData)].min()
        maxdata = ncData[~np.isnan(ncData)].max()
        out_data = []
        #        print(values)
        for k in range(len(values)):
            for j in range(len(values[k])):
                if values[k][j] >= mindata and values[k][j] <= maxdata:
                    out_data.append(values[k][j])
        #        values=np.array(values)
        #        maxs=np.max(values)
        #        means=np.mean(values)
        #        mins=np.min(values)
        values = np.array(out_data)
        maxs = values.max()
        means = values.mean()
        mins = values.min()
        std = values.std()
        xyValues.append([name, maxs, means, mins,std])

    del feature, shpdata

    return xyValues


"""
Calculate the max, min and mean of the raster(image) data according to the polygon, and save to .csv file
inuputShapefile:path of shapefile
inputRasterDir：path of daily environmental factor data, with .nc4 format
(lonResolution,latResolution)：resollution of input environmental factor data
attributeField：name of the input environmental factor, e.g. daily_QV2M
nLon：longititude of environmental factor data
nLat：latitude of environmental factor data
keyField：attribute of Polygon data (shapefile)
dataTypes：used for the names of output results
saveCSVDir：saving path
Dimesion：dimension of input environmental factor data
tempTif：path of temp data
@author: LW, QL
"""

def exportStatisticToCSV(inuputShapefile, inputRasterDir, lonResolution, latResolution, attributeField, nLon, nLat,
                         keyField, dataTypes, saveCSVDir, Dimesion, tempTif):
    all_nc_files, n_files = get_files(inputRasterDir, '*.nc4')
    all_nc_files = np.sort(all_nc_files)
    for i in range(n_files):
        savedir = saveCSVDir
        the_filename = all_nc_files[i]
        date = the_filename.split('.nc4')[0][-8:]
        filename = inputRasterDir + the_filename
        data = extract_h4_by_name(filename, attributeField)
        if Dimesion == '3D':
            data = np.squeeze(data)
        #        data=np.array(data)
        #        data=dat.transpose((1,0))
        lats = extract_h4_by_name(filename, nLat)
        lons = extract_h4_by_name(filename, nLon)
        x_min = lons.min()
        x_max = lons.max()
        y_min = lats.min()
        y_max = lats.max()

        N_Lat = len(lats)
        N_Lon = len(lons)

        spei_ds = gdal.GetDriverByName('Gtiff').Create(tempTif, N_Lon, N_Lat, 1, gdal.GDT_Float32)
        # 3.4 set the range of image visualization
        geotransform = (x_min, lonResolution, 0, y_min, 0, latResolution)
        spei_ds.SetGeoTransform(geotransform)

        # coordinates information
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        spei_ds.SetProjection(srs.ExportToWkt())

        # 3.6 output data
        spei_ds.GetRasterBand(1).WriteArray(data)  #
        spei_ds.FlushCache()  #
        spei_ds = None  #

        rasterdata = rasterio.open(tempTif)
        #        Transform = rasterdata._transform

        # search all the vectors, retuen 2D array and save to .csv file
        csvfile = []
        csvfile = readShapefileFrom2DNCFile(inuputShapefile, rasterdata, lonResolution, latResolution, data, keyField)

        csvfile = np.array(csvfile)
        savedir = savedir + date + '_' + dataTypes + '.csv'
        print(' CSVFile is saving  ', savedir)
        np.savetxt(savedir, csvfile, delimiter=',', fmt='%s')

        del rasterdata


if __name__ == "__main__":
    period = 'history'
    inputRasterDir = r'/Users/sivanuhappy/Documents/George Mason University/COVID-19 Research Project/Daily mean/'+period+'_daily/'
    globalShapeFile = r'/Users/sivanuhappy/Documents/George Mason University/COVID-19 Research Project/Daily mean/global_basemap/global_basemap.shp'
    results = r'/Users/sivanuhappy/Documents/George Mason University/COVID-19 Research Project/Daily mean/Results'
    tempDir = r'/Users/sivanuhappy/Documents/George Mason University/COVID-19 Research Project/Daily mean/Results/temp.tif'
    rastername = 'Tro_NO2'
    exportStatisticToCSV(globalShapeFile, inputRasterDir, 0.25, 0.25, rastername, 'nlon', 'nlat', 'GID_2', rastername,
                         savedirs, '2D', tempDir)

    print('Finish')
