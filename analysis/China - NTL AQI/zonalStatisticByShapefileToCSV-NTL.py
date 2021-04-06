# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 15:42:38 2020

@author: LW
"""


import glob,os,sys
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal,gdal_array,osr, ogr
import rasterio 
from rasterio.mask import mask
from geopandas import *;
import geopandas as geopd
from shapely.geometry import Polygon,Point;

def get_files(dir,ext):
    allfiles=[]
    os.chdir(dir)
    for file in glob.glob(ext):
        allfiles.append(file)
    return allfiles, len(allfiles)

def extract_h4_by_name(filename,dsname):
    h4_data = Dataset(filename)
    print(h4_data)
    ds = np.array(h4_data[dsname][:])
    h4_data.close()
    return ds
"""
遍历所有的矢量要素进行掩膜，返回统计二维数组
"""
def readShapefileFrom2DNCFile(shpdir,datas,xReSolution,yReSolution,ncData,keyFieldName):

    shpdata = geopd.GeoDataFrame.from_file(shpdir)
    #获取要素及要素地理位置
    xyValues = [[keyFieldName,'Max','Mean','Min']]
    #投影变换，使矢量数据与栅格数据投影参数一致
    shpdata=shpdata.to_crs(datas.crs)
    for i in range(0, len(shpdata)):
        # 获取矢量数据的features
        geo = shpdata.geometry[i]
        feature = [geo.__geo_interface__]
        itemindex = np.argwhere(shpdata.columns == keyFieldName)[0][0]
        name=shpdata.iloc[i][itemindex]
#        print(len(shpdata.columns),str(shpdata.columns[itemindex]),name,itemindex)
        #通过feature裁剪栅格影像rasterio.mask.mask 
        out_image, out_transform = mask(datas, feature, all_touched=True, crop=True, nodata=datas.nodata)
        # 掩模得到相交区域
#        out_image, out_transform = mask(datas, geometry, all_touched=True, crop=True, nodata=datas.nodata)
        values = out_image.tolist()
        values = values[0]
        
        #去除异常值
        mindata=ncData[~np.isnan(ncData)].min()
        maxdata=ncData[~np.isnan(ncData)].max()
        out_data=[]
#        print(values)
        for k in range(len(values)):
            for j in range(len(values[k])):
                if values[k][j] >=mindata and values[k][j]<=maxdata:
                    out_data.append(values[k][j])
#        values=np.array(values)
#        maxs=np.max(values)
#        means=np.mean(values)
#        mins=np.min(values)
        values=np.array(out_data)
        maxs=values.max()
        means=values.mean()
        mins=values.min()
        
        xyValues.append([name,maxs,means,mins])
           
    del feature,shpdata
    
    return xyValues
def readShapefileFrom2DTIFFile(shpdir,datas,keyFieldName):

    shpdata = geopd.GeoDataFrame.from_file(shpdir)
    #获取要素及要素地理位置
    xyValues = [[keyFieldName,'Max','Mean','Min']]
    #投影变换，使矢量数据与栅格数据投影参数一致
#    shpdata=shpdata.to_crs(datas.crs)
    for i in range(0, len(shpdata)):
        # 获取矢量数据的features
        geo = shpdata.geometry[i]
        feature = [geo.__geo_interface__]
        itemindex = np.argwhere(shpdata.columns == keyFieldName)[0][0]
        name=shpdata.iloc[i][itemindex]
#        print(len(shpdata.columns),str(shpdata.columns[itemindex]),name,itemindex)
        #通过feature裁剪栅格影像rasterio.mask.mask 
        out_image, out_transform = mask(datas, feature, all_touched=True, crop=True, nodata=datas.nodata)
        # 掩模得到相交区域
#        out_image, out_transform = mask(datas, geometry, all_touched=True, crop=True, nodata=datas.nodata)
        values = out_image.tolist()
        values = values[0]
        
        #去除异常值
        bounds=datas.read(1)
        mindata=bounds[~np.isnan(bounds)].min()
        maxdata=bounds[~np.isnan(bounds)].max()
        out_dataMean=[]
        out_data=[]
        for k in range(len(values)):
            for j in range(len(values[k])):
                if values[k][j] >=mindata and values[k][j]<=maxdata:
                    out_data.append(values[k][j])
                    if values[k][j]>=5:
                        out_dataMean.append(values[k][j])

        values=np.array(out_data)
        out_dataMean=np.array(out_dataMean)
        maxs=values.max()
        means=out_dataMean.mean()
        mins=values.min()
        
        xyValues.append([name,maxs,means,mins])
           
    del feature,shpdata
    
    return xyValues
"""
根据Polygon去统计范围内栅格值的最大最小和均值，并保存到CSV文件中
inuputShapefile:输入的shapefile文件
inputRasterDir：输入的NC栅格文件路径
(lonResolution,latResolution)：输入NC栅格文件分辨率
attributeField：输入NC栅格文件的属性表数组名
nLon：输入NC栅格文件的经度数组名
nLat：输入NC栅格文件的纬度数组名
keyField：Polygon属性表中统计字段名
dataTypes：统计的NC栅格文件对象名
saveCSVDir：保存路径
Dimesion：输入的NC栅格文件维度
tempTif：保存临时文件路径
@author: LW
"""
def exportStatisticToCSV(inuputShapefile,inputRasterDir,lonResolution,latResolution,attributeField,nLon,nLat,keyField,dataTypes,saveCSVDir,Dimesion,tempTif):
    all_nc_files, n_files = get_files(inputRasterDir, '*.nc4')
    all_nc_files = np.sort(all_nc_files)
    for i in range(n_files):
        savedir=saveCSVDir
        the_filename = all_nc_files[i]
        date = the_filename.split('.nc4')[0][-8:]
        filename = inputRasterDir+the_filename
        data = extract_h4_by_name(filename, attributeField)
        if Dimesion=='3D':
            data = np.squeeze(data)
#        data=np.array(data)
#        data=dat.transpose((1,0))
        lats=extract_h4_by_name(filename,nLat)
        lons=extract_h4_by_name(filename,nLon)
        
        lonlat=[]
        for i in range(len(lons)):
            lonlat.append([lons[i],lats[i]])
        
        x_min=lons.min()
        x_max=lons.max()
        y_min=lats.min()
        y_max=lats.max()
        
        N_Lat = len(lats) 
        N_Lon = len(lons)
        
        spei_ds = gdal.GetDriverByName('Gtiff').Create(tempTif,N_Lon,N_Lat,1,gdal.GDT_Float32)
        # 3.4 设置影像的显示范围
        geotransform = (x_min,lonResolution, 0, y_min, 0, latResolution)
        spei_ds.SetGeoTransform(geotransform)
        
        #地理坐标系统信息
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        spei_ds.SetProjection(srs.ExportToWkt())
        
        # 3.6 数据写出
        spei_ds.GetRasterBand(1).WriteArray(data) # 将数据写入内存，此时没有写入硬盘
        spei_ds.FlushCache() # 将数据写入硬盘
        spei_ds = None # 关闭spei_ds指针，注意必须关闭
        
        rasterdata = rasterio.open(tempTif)
#        Transform = rasterdata._transform  # 得到影像六参数
        
        #遍历所有的矢量要素，返回二维数组，保存到csv文件中
        csvfile=[]
        csvfile=readShapefileFrom2DNCFile(inuputShapefile,rasterdata,lonResolution,latResolution,data,keyField)
        
        csvfile=np.array(csvfile)
        savedir=savedir+date+'_'+dataTypes+'.csv'
        print(' CSVFile is saving  ',savedir)
        np.savetxt(savedir,csvfile,delimiter = ',', fmt='%s')
        
        del rasterdata

def exportTIFStatisticToCSV(inuputShapefile,inputRasterDir,lonResolution,latResolution,attributeField,nLon,nLat,keyField,dataTypes,saveCSVDir,Dimesion,tempTif):
    all_nc_files, n_files = get_files(inputRasterDir, '*.tif')
    all_nc_files = np.sort(all_nc_files)
    for i in range(n_files):
        savedir=saveCSVDir
        inputdir=inputRasterDir+'/'+all_nc_files[i]
        rasterdata = rasterio.open(inputdir)
#        Transform = rasterdata._transform  # 得到影像六参数
#        print(rasterdata.read(1).max(),rasterdata.read(1).min())
        #遍历所有的矢量要素，返回二维数组，保存到csv文件中
        csvfile=[]
        csvfile=readShapefileFrom2DTIFFile(inuputShapefile,rasterdata,keyField)

        csvfile=np.array(csvfile)
        savedir=savedir+'_'+dataTypes+'.csv'
        print(' CSVFile is saving  ',savedir)
        np.savetxt(savedir,csvfile,delimiter = ',', fmt='%s')
        
        #进行统计范围个数
        band=rasterdata.read(1)
        
        band[band>40]=-1
        band[band>20]=-2
        band[band>5]=-3
        band[band>0]=-4
        range520=0
        range2040=0
        range40=0
#        key = np.unique(band)
        key=[-1,-2,-3,-4]
        print(key)
        result = [['Range','count']]
        
        for k in key:
            mask = (band == k)
            band_new = band[mask]
            v = band_new.size
            if k==-3:
                range520=range520+v
            if k==-2:
                range2040=range2040+v
            if k==-1:
                range40=range40+v    
            
        result.append(['5~20',range520])
        result.append(['20~40',range2040])
        result.append(['>40',range40])
        savedir=saveCSVDir+'_'+dataTypes+'Statistic.csv'
        print(' CSVFile is saving  ',savedir)
        np.savetxt(savedir,result,delimiter = ',', fmt='%s')
        del rasterdata

if __name__ == "__main__": 
    
    # 定义数据路径

    inputPrecipRasterDir=r'F:/GMU-COVID-19/experiment/experiment/precipitation-new/'
    inputRHRasterDir=r'F:/GMU-COVID-19/experiment/experiment/humidity/'
    
    
    inputTempRasterDir=r'F:/GMU-COVID-19/qliu/'
    
    admin2_USA=r'F:\GMU-COVID-19\qliu\CHN_shp-province\CHN_shp\gadm36_CHN_1.shp'

    savedirs=r'F:/GMU-COVID-19/CSVResult/'
    
    tempDir=r'F:/GMU-COVID-19/CSVResult/temp.tif'
    
    #开始执行统计，并存储CSV
#    dataType='Humidity'
#    exportStatisticToCSV(inputRHRasterDir,0.625,0.5,'daily_QV2M','nlon','nlat','GID_2',dataType,savedirs,'2D',tempDir)
#    
#    dataType='Precipitation'
#    exportStatisticToCSV(inputPrecipRasterDir,0.1,0.1,'daily_precipitation','nlon','nlat','GID_2',dataType,savedirs,'2D',tempDir)
#    
    dataType='Radiance_difference'
    exportTIFStatisticToCSV(admin2_USA,inputTempRasterDir,0.625,0.5,'radiance difference','nlon','nlat','GID_1',dataType,savedirs,'2D',tempDir)
    
    
 
    print('Finish')   