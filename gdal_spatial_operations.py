import os
import sys
from osgeo import gdal, ogr, osr
import re
import math
import pandas as pd

from gdal_projection import *

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()




def vector_to_geotiff(inVector, AttributeField, outRaster, pixel_size = None, NoData = -9999, DType = 'GDT_byte', snap = False, snap_raster = None, shape_type = None):
    '''
    Function Def: vector_to_geotiff(inVector, AttributeField, outRaster, pixel_size, NoData = -9999, DType = 'GDT_byte', snap = False, shape_type = None)

    inVector        ~ Input vector file (shapefile), it can be either point line or polygon [REQUIRED]
    AttributeField  ~ The exact name of the column that the rasterization will be based of (in the attibute table) NB: This needs to be a NUMERICAL column at this point  [REQUIRED]
    outRaster       ~ Name of the output raster (path and name, follwed by .tif) [REQUIRED]
    pixel_size      ~ Pixel size of the output raster in meters (for both geographic and projected coordinates), if snape is true, the value that you have entered will be overwritten [IF NOT SNAP REQUIRED]
    NoData          ~ NoData value of the output raster [OPTIONAL]
    DType           ~ Data type out the output raster, Choises are exclusively: GDT_byte; GDT_UInt16; GDT_Int16; GDT_Float32 [OPTIONAL]
    snap            ~ Whether to snap the output to the snap_raster or not, exclusively True or False [OPTIONAL]
    snap_raster     ~ The raster to which the output will be snapped, IF snap is True, this field is [IF SNAP REQUIRED]
    shape_type      ~ exclusively either "point", "line" or "polygon" of inVector, this is only needed if snape is true [IF SNAP REQUIRED]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns either: ~ if snap is False, a raster that represents the vector file in the vector's coordinate system
                    ~ if snape is True, a snapped raster that represents the vector file in the snap raster's coordinate system
    '''       
    try:
        
        shp_temp_status = 'No'
        if snap == True:
            if snap_raster == None:
                raise Exception('Error: Vector to geotiff: Snap is indicated as True without a snap raster')
            else:
                SNR = gdal.Open(snap_raster)        # Getting snap raster details
                gt_snap = SNR.GetGeoTransform()     # ulx, xres, xskew, uly, yskew, yres
                xRS = x_arr_size = SNR.RasterXSize
                yRS = y_arr_size = SNR.RasterYSize

                x_pix_size = gt_snap[1]             # Setting up pixel size and extents
                y_pix_size = gt_snap[5]
                x_min = gt_snap[0]
                y_max = gt_snap[3]
                x_max = x_min + (x_pix_size*xRS)
                y_min = y_max + (y_pix_size*yRS)
                
                vector_coord = get_vector_epsg(inVector)
                vector_coord = vector_coord[0]
                raster_coord = get_raster_epsg(snap_raster)
                raster_coord = raster_coord[0]
                
                if vector_coord != raster_coord:
                    temp_vector = inVector[0:len(inVector)-4] + '_temp.shp'
                    project_shp(inVector, temp_vector, raster_coord, shape_type)
                    inVector = temp_vector
                    shp_temp_status = 'Yes'
        else:
            if pixel_size == None:
                raise Exception('Error: Vector to geotiff: Snap is indicated as False and no output pixel size is defined')
 
        
        vector_DS = ogr.Open(inVector)                          # [IN] Opening the input vector filer [Driver]->[Datasource]
        vector_layer = vector_DS.GetLayer()                     # [IN] Creating layer object from the datasource object [Datasource]->[Layer]

        CStype = wkt2epsg(get_shapefile_wkt(inVector), cstypeStat = True)   # [IN] Checks the coordinate system type of the input shapefile
        # If the inputfile is in geographic coordinates, the pixel cell size is determined by converting it to a projected coordinate system
        if CStype == 'GEOGCS' and snap == False:
            #________________________________________________________looks for the UTM zone EPSG code (pretty schweet)
            x_min, x_max, y_min, y_max = vector_layer.GetExtent()   # [IN] Get the extent of the vector file (Extent coordinates will be in the format of the input vector files projection)[Layer]
            ePsG_list = get_vector_epsg(inVector)
            inEpsg = ePsG_list[0]
            inEpsg_num = int(inEpsg[5:len(inEpsg)])
            inSpatial_ref = osr.SpatialReference()
            inSpatial_ref.ImportFromEPSG(inEpsg_num)
            utm_band = str((math.floor((x_min + 180)/6)%60)+1)
            if  len(utm_band) == 1:
                utm_band = '0'+ utm_band
            if y_max >= 0:
                epsg_code = '326' + utm_band
            else:
                epsg_code = '327' + utm_band
            epsg_code = int(epsg_code)
            outSpatial_ref = osr.SpatialReference()
            outSpatial_ref.ImportFromEPSG(epsg_code)
            Geo_to_Proj = osr.CoordinateTransformation(inSpatial_ref, outSpatial_ref)        
            (ulx, uly, ulz) = Geo_to_Proj.TransformPoint(x_min, y_max)
            (lrx, lry, lrz) = Geo_to_Proj.TransformPoint(x_max, y_min)
            num_of_pixels = int(round((lrx- ulx)/pixel_size, 0))
            pixel_size = abs((x_max - x_min)/num_of_pixels)
            x_pix_size = pixel_size                             
            y_pix_size = -pixel_size
            x_arr_size = abs(int(round((x_max - x_min)/x_pix_size,0))) + 2          # Computing the x size of the output Gtiff array
            y_arr_size = abs(int(round((y_max - y_min)/abs(y_pix_size),0))) + 2     # Computing the y size of the output Gtiff array
            
            
        elif CStype == 'PROJCS' and snap == False:
            x_pix_size = pixel_size
            y_pix_size = -pixel_size
            x_min, x_max, y_min, y_max = vector_layer.GetExtent()   # [IN] Get the extent of the vector file (Extent coordinates will be in the format of the input vector files projection)[Layer]
            x_arr_size = abs(int(round((x_max - x_min)/x_pix_size,0))) + 2          # Computing the x size of the output Gtiff array
            y_arr_size = abs(int(round((y_max - y_min)/abs(y_pix_size),0))) + 2     # Computing the y size of the output Gtiff array
            

        if DType == 'GDT_byte':
            output_DS = gdal.GetDriverByName('MEM').Create('memory_raster', x_arr_size, y_arr_size,1 ,gdal.GDT_Byte)    #[OUT] Creates an empty Gtiff raster with the x and y sizes calculated in memory [Driver]->[Dataset]
        elif DType == 'GDT_UInt16':
            output_DS = gdal.GetDriverByName('MEM').Create('memory_raster', x_arr_size, y_arr_size,1 ,gdal.GDT_UInt16)  #[OUT] Creates an empty Gtiff raster with the x and y sizes calculated in memory [Driver]->[Dataset]
        elif DType == 'GDT_Int16':
            output_DS = gdal.GetDriverByName('MEM').Create('memory_raster', x_arr_size, y_arr_size,1 ,gdal.GDT_Int16)   #[OUT] Creates an empty Gtiff raster with the x and y sizes calculated in memory [Driver]->[Dataset]
        elif DType == 'GDT_Float32':
            output_DS = gdal.GetDriverByName('MEM').Create('memory_raster', x_arr_size, y_arr_size,1 ,gdal.GDT_Float32) #[OUT] Creates an empty Gtiff raster with the x and y sizes calculated in memory [Driver]->[Dataset]

        output_DS.SetProjection(get_shapefile_wkt(inVector))
        output_DS.SetGeoTransform((x_min, x_pix_size, 0, y_max, 0, y_pix_size))  # [OUT] Setting geotransformation of the empty Gtiff array [Dataset] # ulx, xres, xskew, uly, yskew, yres
        

        output_Band = output_DS.GetRasterBand(1)    # [OUT] Creates a band object from the empty raster in memory [Dataset]->[Band]
        output_Band.SetNoDataValue(NoData)          # [OUT] Sets the NoData value of the band

        nodata_array = np.full((y_arr_size, x_arr_size), NoData)    # NoData HACK
        output_Band.WriteArray(nodata_array)                        # NoData HACK
        
    except Exception as err1:
        raise Exception('Error in processing data for Feature to Raster conversion:', err1)

    try:
        gdal.RasterizeLayer(output_DS, [1], vector_layer, options = ['ATTRIBUTE=%s' % AttributeField]) # [OUT] [Function] rasterizes the vector file into the created output file using the specified attribute field
    except Exception as err2:
        raise Exception('Error in Atempting to rasterize feature layer: ', err2)
        
    gdal.GetDriverByName('GTiff').CreateCopy(outRaster, output_DS)

    if shp_temp_status == 'Yes':
        vector_DS = None
        driver = ogr.GetDriverByName('ESRI Shapefile')
        driver.DeleteDataSource(inVector)            
         


def pandas_df_to_point(DB_query, z_val, outShape, limit_extent = None):
    '''
    pandas_df_to_point(DB_query, z_val, outShape, limit_extent = None):
    STILL TO DO
    '''
    
    if limit_extent != None:                                    #### FOR EXTENT LIMITER ####                 
        if os.path.isfile(limit_extent):
            coord_pair = get_raster_epsg(limit_extent)
            extent_coord = coord_pair[0]
            if extent_coord != 'EPSG:4326':
                project_raster(limit_extent, limit_extent[0:len(limit_extent)-4] + '_Temp_extent.tif', 'EPSG:4326')
                extR = limit_extent[0:len(limit_extent)-4] + '_Temp_extent.tif'
                extent_rst = gdal.Open(extR)
                extent_gt = extent_rst.GetGeoTransform()     # Geotransformation # ulx[0], xres[1], xskew[2], uly[3], yskew[4], yres[5]
                xRS =  extent_rst.RasterXSize                   # Number of columns
                yRS =  extent_rst.RasterYSize                   # Number of rows
                extent_rst = None
                driver_del = gdal.GetDriverByName('GTiff')
                driver_del.Delete(limit_extent[0:len(limit_extent)-4] + '_Temp_extent.tif')
                
            else:
                extR = limit_extent
                extent_rst = gdal.Open(extR)
                extent_gt = extent_rst.GetGeoTransform()     # Geotransformation # ulx[0], xres[1], xskew[2], uly[3], yskew[4], yres[5]
                xRS =  extent_rst.RasterXSize                   # Number of columns
                yRS =  extent_rst.RasterYSize                   # Number of rows
                extent_rst = None
                
            Xul = extent_gt[0]
            Xlr = Xul + (xRS * extent_gt[1])
            Yul = extent_gt[3]
            Ylr = Yul + (yRS * extent_gt[5])                # extent_gt[5] will be minus

        else:
            raise Exception ('ERROR in converting data base to points: Extent limit file specified does not exist')


    driver = ogr.GetDriverByName('ESRI Shapefile')

    if os.path.exists(outShape):                            # [OUT] Checks whether the output shapefile already exists, if it does, it gets deleted
        driver.DeleteDataSource(outShape)    
    outDataSet = driver.CreateDataSource(outShape)

    spatialRef = osr.SpatialReference()                             # [PROJ][OSR] Since an esri file is created the .prj file still has to be written in the suitable format. This lines creates another spatial reference object
    spatialRef.ImportFromEPSG(4326)
    
    outLayer = outDataSet.CreateLayer(z_val, spatialRef, ogr.wkbPoint)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    attribute = ogr.FieldDefn(z_val[0:10], ogr.OFTReal)
    outLayer.CreateField(attribute)
     
    for index, row in DB_query.iterrows():

        if limit_extent != None:

            X_p = float(row[0])
            Y_p = float(row[1])

            if X_p > Xul and X_p < Xlr and Y_p < Yul and Y_p > Ylr:         # extent query to see whether the point falls within the specified extent
                feature = ogr.Feature(outLayer.GetLayerDefn())
                feature.SetField(z_val[0:10], row[2])
                wkt = "POINT(%f %f)" % (float(row[0]), float(row[1]))
                point = ogr.CreateGeometryFromWkt(wkt)
                feature.SetGeometry(point)
                outLayer.CreateFeature(feature)
                feature = None
            else:
                DB_query.drop(index, inplace=True)

        else:
            feature = ogr.Feature(outLayer.GetLayerDefn())
            feature.SetField(z_val[0:10], row[2])
            wkt = "POINT(%f %f)" % (float(row[0]), float(row[1]))
            point = ogr.CreateGeometryFromWkt(wkt)
            feature.SetGeometry(point)
            outLayer.CreateFeature(feature)
            feature = None

    outDataSet = None

    return DB_query



         

'''
out = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/'

point = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/TMP.shp'
point_prj = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/TMP_PRJ.shp'
DEM_wgs = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_wgs.tif'
DEM_UTM = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_prj_UTM_Zone_34S.tif'

print('vector to raster version 1')
vector_to_geotiff(point, 'Temp', out + 'temp_ras_geo_VR1.tif', 30, DType = 'GDT_Float32')
print('vector to raster version 2')
vector_to_geotiff(point_prj, 'Temp', out + 'temp_ras_geo_VR2.tif', 30, DType = 'GDT_Float32')
print('vector to raster version 3')
vector_to_geotiff(point, 'Temp', out + 'temp_ras_geo_VR3.tif', DType = 'GDT_Float32', snap = True, snap_raster = DEM_wgs, shape_type = 'point')
print('vector to raster version 4')
vector_to_geotiff(point, 'Temp', out + 'temp_ras_geo_VR4.tif', DType = 'GDT_Float32', snap = True, snap_raster = DEM_UTM, shape_type = 'point')
print("DONE")
'''
#____________________________
'''
polygon = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/square.shp'
polygonOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/squareRAS_GEOG.tif'

point = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/TMP_PRJ.shp'
pointOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/TMP_RAS.tif'

line = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/test_linePRJ.shp'
lineOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/test_line_RAS.tif'


vector_to_geotiff(polygon, 'rasterize', polygonOUT, 30)
#vector_to_geotiff(point, 'Temp', pointOUT, 30, DType = 'GDT_Float32')
#vector_to_geotiff(line, 'rasterize', lineOUT, 30)
print('DONE')
#vector_to_geotiff('vector', 'Attribur', 30, 34234)
'''
