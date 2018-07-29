import os
import sys
from osgeo import gdal, ogr, osr
import re
import math
import numpy as np

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()



def wkt2epsg(wkt, forceProj4 = False, cstypeStat = False):
    '''
    Function Def: wkt2epsg(wkt, forceProj4 = False, cstypeStat = False)

    Adapted from source:
    source: https://gis.stackexchange.com/questions/20298/is-it-possible-to-get-the-epsg-value-from-an-osr-spatialreference-class-using-th
    
    wkt         ~  Porjection information in the wkt format [REQUIRED]
    forceProj4  ~  This is if there is a problem and you want to force a Proj4 text as a return value [OPTIONAL]
    cstypeStat  ~  If this is specified as True, the function will only return if the projection is a Geographic or a projected coordinate system. [OPTIONAL]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    
    returns either:    ~ list containing epsg code, if found, and a confidence ranking (1 to 6, with 1 being the highes confidence) .e.g ['EPSG:32734', 1]
                       ~ None, if the EPSG code could not be found.
                       ~ cstype is True, the function will only return 'GEOGCS' or 'PROJCS'                 
    '''
                       
    rel_dir = os.path.dirname(__file__) # path of the epsg is always relative to this file
    rel_dir = rel_dir + '/'
    epsg = rel_dir + 'epsg'
    esri = rel_dir + 'esri'
    code= None
    
    p_in = osr.SpatialReference()
    s = p_in.ImportFromWkt(wkt)

    if s == 5:  # invalid WKT
        return None
    if p_in.IsLocal() == 1:             # this is a local definition
        return p_in.ExportToWkt()
    if p_in.IsGeographic() == 1:        # This is a geographic srs
        cstype = 'GEOGCS'
    else:
        cstype = 'PROJCS'               # This is a projected srs

    if cstypeStat == True:
        return cstype
    
    an = p_in.GetAuthorityName(cstype)
    ac = p_in.GetAuthorityCode(cstype)

    if an is not None and ac is not None:   # return EPSG code
        return ['%s:%s' % (p_in.GetAuthorityName(cstype), p_in.GetAuthorityCode(cstype)),1]
    else: # try brute force appproach by grokking proj epsg definition file
        p_out = p_in.ExportToProj4()
        
        if p_out:
            if forceProj4 is True:
                return [p_out,6]

            f = open(epsg)
            for line in f:       
                if line.find(p_out) != -1:
                    m = re.search('<(\\d+)>', line)
                    if m:
                        code = m.group(1)
                        break                     
            if code:
                return ['EPSG:%s' % code,2]
            else:               # did not find the spatial reference in the epsg file, now systematically searching the esri file which has a slightly different format
                p_out_ellps = None
                p_out_datum = None
                   
                if p_out.find('+ellps=') != -1:
                    elp_res = re.search(r'\+ellps=(\S+)', p_out)
                    p_out_ellps = elp_res.group()
                    p_out_ellps = p_out_ellps.strip()

                if p_out.find('+datum=') != -1:
                    dt_res = re.search(r'\+datum=(\S+)', p_out)
                    p_out_datum = dt_res.group()
                    p_out_datum = p_out_datum.strip()
                    
                if p_out_ellps == None and p_out_datum == None:
                    if p_out.find('+no_defs') != -1:
                        nodef_1 = re.search(r'(.+)(?=\+no_def)', p_out)
                        new_p_out = nodef_1.group(1)
                        new_p_out = new_p_out.strip()
                    else:
                        nodef_1 = re.search(r'(.+)(?=no_defs)', p_out)
                        new_p_out = nodef_1.group(1)
                        new_p_out = new_p_out.strip()
                        
                elif p_out_ellps == None and p_out_datum != None:
                    rslt = re.search(r'(.+)(?=\+datum=)', p_out)
                    new_p_out = rslt.group(1)
                    new_p_out = new_p_out.strip()
                else:
                    rslt = re.search(r'(.+)(?=\+ellps=)', p_out)
                    new_p_out = rslt.group(1)
                    new_p_out = new_p_out.strip()
                
                f.close()
                f = open(esri)
                for line in f:

                    if  line[0] != '#' and (line[0:8] == '<102022>' or line[0:6] == '<4326>' or line[0:7] == '<32734>'): 

                        line_ellps = None
                        line_datum = None
                        spatial_code = 'None'   
                        if line.find('+ellps=') != -1:
                            elp_res = re.search(r'\+ellps=(\S+)', line)
                            line_ellps = elp_res.group()

                        if line.find('+datum=') != -1:
                            dt_res = re.search(r'\+datum=(\S+)', line)
                            line_datum = dt_res.group()

                        if line_ellps == None and line_datum == None:
                            if line.find('+no_defs') != -1:
                                nodef_1 = re.search(r'<\d+>(.+)(?=\+no_def)', line)
                                new_line = nodef_1.group(1)
                                new_line = new_line.strip()
                            else:
                                nodef_1 = re.search(r'<\d+>(.+)(?=no_defs)', line)
                                new_line = nodef_1.group(1)
                                new_line = new_line.strip()
                                
                        elif line_ellps == None and line_datum != None:
                            rslt = re.search(r'<\d+>(.+)(?=\+datum=)', line)
                            new_line = rslt.group(1)
                            new_line = new_line.strip()
                        else:
                            rslt = re.search(r'<\d+>(.+)(?=\+ellps=)', line)
                            new_line = rslt.group(1)
                            new_line = new_line.strip()

                        if len(new_line) == 0:      # This is not 100% water tight becuase the know EPSG code are usually catched in either the first of second phases, however more logical checks can be placed here
                            new_line = 'NoThInG'
                        elif new_line == '+proj=longlat':
                            if line_ellps != p_out_ellps and line_datum != p_out_datum:
                                new_line = 'NoThInG'
                                                    
                        if new_line == new_p_out:
                            if (p_out_datum == line_datum and p_out_datum != None and line_datum != None) and \
                            (p_out_ellps == line_ellps and p_out_ellps != None and line_ellps != None):
                                m = re.search('<(\\d+)>', line)
                                if m:
                                    code = m.group(1)
                                    spatial_code = 'D&E'
                                    break
                            elif (p_out_datum == line_datum and p_out_datum != None and line_datum != None) or \
                            (p_out_ellps == line_ellps and p_out_ellps != None and line_ellps != None):
                                m = re.search('<(\\d+)>', line)
                                if m:
                                    code = m.group(1)
                                    spatial_code = 'D|E'
                                    break
                            else:
                                m = re.search('<(\\d+)>', line)
                                if m:
                                    code = m.group(1)
                                    spatial_code = '!D!E'
                                    break
                            
                if code:
                    if spatial_code == 'D&E':                       
                        return ['EPSG:%s' %code, 3]
                    elif spatial_code == 'D|E':
                        return ['EPSG:%s' %code, 4]
                    elif spatial_code == '!D!E':
                        return ['EPSG:%s' %code, 5]
                else:
                    return None
        else:
            return None
                
def get_raster_wkt(input_raster):
    '''
    Function Def: get_raster_wkt(input_raster)
    
    input_raster ~ the path and name of the input raster image (e.g. .tif) [REQUIRED]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns     ~ The function returns the spatial reference object of input_raster
                  in an well known text (WKT) format
    '''
    input_layer = gdal.Open(input_raster)
    prj_ras = input_layer.GetProjection()
    srs_ras = osr.SpatialReference(wkt=prj_ras)
    return srs_ras.ExportToWkt()

def get_shapefile_wkt(input_shape):
    '''
    Function Def: get_shapefile_wkt(input_shape)
    
    input_shape ~ the ath and name of the input shapefile(.shp) [REQUIRED]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns     ~ The function returns the spatial reference object of input_shape
                  in an well known text (WKT) format
    '''
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataset = driver.Open(input_shape)
    layer = dataset.GetLayer()
    spatialRef = layer.GetSpatialRef()
    return spatialRef.ExportToWkt()

def project_shp(inShape, outShape, outEpsg, shape_type): 

    '''
    Function Def: project_shp(inShape, outShape, outEpsg, shape_type): 
    
    inShape      ~ the path and name of the input shapefile (.shp) [REQUIRED]
    outShape     ~ the path and name of your desired output projected shapefile (.shp) [REQUIRED]
    outEpsg_num  ~ the output coordinate system in EPSG format e.g. "EPSG:4326" [REQUIRED]
    shape_type   ~ exclusively either "point", "line" or "polygon" [REQUIRED]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns      ~ This function does not return anything but outputs a projected vector file at outShape
    '''

    ePsG_list = wkt2epsg(get_shapefile_wkt(inShape))    # [PROJ] Gets the EPSG number

    if ePsG_list == None:
        raise Exception('PROJECTION ERROR: The input coordinate system (EPSG code) could not be defined. The file most likely contains a custom projection made in a geospatial program without an known EPSG code\n'  + \
                        'Please re-project the file: ' + inShape + ' to a defined coordinate system then try again')
                                                        
    inEpsg = ePsG_list[0]                               # [PROJ][OSR] Extracting the EPSG numbers 
    Epsg_num = int(inEpsg[5:len(inEpsg)])
    outEpsg_num = int(outEpsg[5:len(outEpsg)])

    inSpatialRef = osr.SpatialReference()               # [PROJ][OSR] Creates a spatial reference object for the Input Spatial Reference 
    inSpatialRef.ImportFromEPSG(Epsg_num)               # [PROJ][OSR] Populates the spatial refernce object with the input EPSG number 
        
    outSpatialRef = osr.SpatialReference()              # [PROJ][OSR] Creates a spatial reference object for the Output Spatial Reference
    outSpatialRef.ImportFromEPSG(outEpsg_num)           # [PROJ][OSR] Populates the spatial reference object with the user defined EPSG number
    
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef) # Creates a coordinate transform object that transforms the input coordinates to the output coordinates
    #_____________________________________________________________________
    
    driver = ogr.GetDriverByName('ESRI Shapefile')      # [IN/OUT]Creating a ESRI Driver object, this Driver will be used throughout [Driver]
    #_____________________________________________________________________

    inDataSet = driver.Open(inShape)                    # [IN] Opening the input vector file with the Dirver object to create a Datasource object  [Driver]->[Datasource]
    inLayer = inDataSet.GetLayer()                      # [IN] Creates a Layer object from the Datasource object [Datasource]->[Layer]
    num_of_in_ftrs = inLayer.GetFeatureCount()          # [IN] Gets the amount of features (lines points or polygons) of the input feature [Layer] ##(Not used in the script, but kept as an example)
    inLayerDefn = inLayer.GetLayerDefn()                # [IN] Gets information about features. Encapsulates the attribute schems of the features of the layer [Layer]->[FeatureDefn]
    attribute_fields = inLayerDefn.GetFieldCount()      # [IN] Number of attributes fields in the input vector [FeatureDefn]

    #_____________________________________________________________________
    
    if os.path.exists(outShape):                            # [OUT] Checks whether the output shapefile already exists, if it does, it gets deleted
        driver.DeleteDataSource(outShape)    
    outDataSet = driver.CreateDataSource(outShape)          # [OUT] Creates a new empty Datasource object by using the Driver object, this will be the output dataset [Driver]->[Datasource] 

    if shape_type == 'polygon' or shape_type == 'Polygon':  # [OUT] Creates a new empty Layer object for the output Datasource object depending on the vector type [Datasource]->[Layer]
        outLayer = outDataSet.CreateLayer('My_layer', geom_type=ogr.wkbMultiPolygon)    

    elif shape_type == 'point'or shape_type == 'Point':
        outLayer = outDataSet.CreateLayer('My_layer', geom_type=ogr.wkbPoint)

    elif shape_type == 'line' or shape_type == 'Line':
        outLayer = outDataSet.CreateLayer('My_layer', geom_type=ogr.wkbLineString)
        
    for i in range(0,attribute_fields):             # For each attirbute field as defined by the input vector field          
        fieldDefn = inLayerDefn.GetFieldDefn(i)     # [IN/OUT] Gets the sepcific information of attribute fiel i [FeatureDefn]
        outLayer.CreateField(fieldDefn)             # [OUT] Adds the field definition to the created output Layer [FeatureDefn]->[Layer]
        
    outLayerDefn = outLayer.GetLayerDefn()          # [OUT] Creates a FeatureDefn object of the outLayer [Layer]->[FeatureDefn]

    #_____________________________________________________________________
    
    inFeature = inLayer.GetNextFeature()        # [IN] Itterates through each feature within the input Layer [Layer]->[Feature]
    while inFeature:
        geom = inFeature.GetGeometryRef()       # [IN] Gets the geometry from the current input feature [Feature]->[Geometry]
        geom.Transform(coordTrans)              # [PROJ] Uses OSR coordinate transform object to convert the geometry object to the new coordinate system [Geometry]

        outFeature = ogr.Feature(outLayerDefn)  # [OUT] Creates an empty feature object and the feature definition of the output layer [Feature]
        outFeature.SetGeometry(geom)            # [OUT] Sets the geometry of the newly created feature with the projected geometry object [Feature]

        for i in range(0, outLayerDefn.GetFieldCount()):    # [OUT] For the length of the output layer definition (would be the same as for the input atrribute fields)
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i)) # Sets the atttribute information in the newlt created feature [Feature]

        try:
            outLayer.CreateFeature(outFeature)      # [OUT] Adds the newly created feature to the created layer [Feature]->[Layer]
        except RuntimeError:
            raise Exception('PROJECTION ERROR: At this time, only supports: wkbMultiPolygon; wkbPoint; wkbLineString')
        outFeature = None                       # [OUT] Cleans the feature so that a new feature can be created again [Feature]
        inFeature = inLayer.GetNextFeature()    # [IN] Moves to the next input feautre (ittirative process) [Feature]

    inDataSet = None        # [IN] At this point the new data set is created and the input dataset can be cleared [Datasource]
    outDataSet = None       # [OUT] At this point the new data set is created and the newly created output dataset can be closed [Datasource]
    #_____________________________________________________________________

    spatialRef = osr.SpatialReference()         # [PROJ][OSR] Since an esri file is created the .prj file still has to be written in the suitable format. This lines creates another spatial reference object
    spatialRef.ImportFromEPSG(outEpsg_num)      # [PROJ][OSR] Populates the spatial refernce object with the input EPSG number 
    spatialRef.MorphToESRI()                    # [PROJ][OSR] Morphs the spatial reference object to a format ESRI software (and most geospatial sofwares) will understand.
    file = open(outShape[0:len(outShape)-4] + '.prj', 'w+') #[PROJ][OSR] Writes the spatial refence object to a .prj file Well Known Text format. 
    file.write(spatialRef.ExportToWkt())
    file.close()

def project_raster(inRas, outRas, outEpsg, outPixelsize=None, NoDataValue = -9999, snap_clip=False, Snap_Raster=None, resample_alg = 'Nearest', data_type = None):   
    '''
    Function Def: project_raster(inRas, outRas, outEpsg, outPixelsize=None, NoDataValue = -9999, snap_clip=False, Snap_Raster=None, resample_alg = 'Nearest')
    
    inRas        ~ the path and name of the input raster (e.g. .tif) [REQUIRED]
    outRas       ~ the path and name of the desired output projected raster [REQUIRED]
    outEpsg      ~ the output coordinate system in EPSG formate.g. 'EPSG:4326', this will be overwritten if snap is True [REQUIRED]
    outPixelsize ~ the pixel size of outRas in meters (the default is its own pixel size), if snap is true, the pixel size will be overwritten with snap raster pixel size [OPTIONAL]
                 ~ if you specify the pixelsize to somehting other than its owns pixel size, the function can serve as a resampling fucntion
    NoDataValue  ~ The Nodata value of the output raster
    snap_clip    ~ weather or not to run the snap raster process, the default is FALSE, the pixel size and extent of the snap raster will overwrite outPixelsize [OPTIONAL]
    resample_alg ~ What algorithm to use for the resampling process, the default is Nearest Neigbour, can be, exclusively: 'Nearest', 'Bilinear', 'Cubic', 'Cubic_spline' [OPTIONAL]

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    returns      ~ This function does not return anything, but outputs a projected raster, either snapped or not at outRas
   
    '''
    
    ePsG_list = get_raster_epsg(inRas)

    if ePsG_list == None:
        raise Exception('PROJECTION ERROR: The input coordinate system (EPSG code) could not be defined. The file most likely contains a custom projection made in a geospatial program without an known EPSG code\n'  + \
                        'Please re-project the file: ' + inRas + ' to a defined coordinate system then try again')
        
    
    epsg_from = ePsG_list[0]
    epsg_from_num = int(epsg_from[5:len(epsg_from)])
    from_coord = osr.SpatialReference()
    from_coord.ImportFromEPSG(epsg_from_num)

    if snap_clip == True:                           ##### STILL SORT OUT DATA TYPE 8BIT 32BIT 
        SnP_epsg = get_raster_epsg(Snap_Raster)
        SnP_epsg = SnP_epsg[0]
        epsg_to_num = int(outEpsg[5:len(SnP_epsg)]) 
        to_coord = osr.SpatialReference()
        to_coord.ImportFromEPSG(epsg_to_num)
        coord_trans = osr.CoordinateTransformation(from_coord, to_coord)

        inData = gdal.Open(inRas)
        gt = inData.GetGeoTransform()
        xRS = inData.RasterXSize
        yRS = inData.RasterYSize

        snapRaster = gdal.Open(Snap_Raster)
        SN_xRS = snapRaster.RasterXSize
        SN_yRS = snapRaster.RasterYSize
        gtSNAP = snapRaster.GetGeoTransform()   # Geotransformation # ulx[0], xres[1], xskew[2], uly[3], yskew[4], yres[5]         
        
        outPixelsize_x = gtSNAP[1]
        outPixelsize_y = gtSNAP[5]
        
        proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", SN_xRS, SN_yRS, 1, gdal.GDT_Float32)
        proj_ras.SetGeoTransform(gtSNAP)
        proj_ras.SetProjection(to_coord.ExportToWkt())
        proj_ras_band = proj_ras.GetRasterBand(1)       
        proj_ras_band.SetNoDataValue(NoDataValue)   

        nodata_array = np.full((SN_yRS, SN_xRS), NoDataValue)   # NoData HACK
        proj_ras_band.WriteArray(nodata_array)                  # NoData HACK
        
        if resample_alg == 'Nearest':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_NearestNeighbour) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Bilinear':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Bilinear)
        elif resample_alg == 'Cubic':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Cubic)
        elif resample_alg == 'Cubic_spline':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_CubicSpline)
            
        
        gdal.GetDriverByName('GTiff').CreateCopy(outRas, proj_ras)
        inData=None
        proj_ras=None
        
    
    else:
        epsg_to_num = int(outEpsg[5:len(outEpsg)])  
        to_coord = osr.SpatialReference()
        to_coord.ImportFromEPSG(epsg_to_num)

        inData = gdal.Open(inRas)
        gt = inData.GetGeoTransform()
        xRS = inData.RasterXSize
        yRS = inData.RasterYSize

        coord_trans = osr.CoordinateTransformation(from_coord, to_coord)

        (ulx, uly, ulz) = coord_trans.TransformPoint(gt[0], gt[3])
        (lrx, lry, lrz) = coord_trans.TransformPoint(gt[0] + gt[1]*xRS, gt[3] + gt[5]*yRS)

        inCStype = wkt2epsg(get_raster_wkt(inRas), cstypeStat = True)
        outCStype = wkt2epsg(to_coord.ExportToWkt(), cstypeStat = True)
        if outPixelsize is None:    
            if inCStype == outCStype:
                outPixelsize = gt[1]
            elif inCStype == 'GEOGCS' and outCStype == 'PROJCS':
                outPixelsize = int((lrx - ulx)/xRS)
            else:
                outPixelsize = (lrx - ulx)/xRS
        else:
            if outCStype == 'GEOGCS':
                if inCStype == 'GEOGCS':

                    ulx_trns = gt[0]
                    uly_trns = gt[3]

                    inSpatial_ref = osr.SpatialReference()
                    inSpatial_ref.ImportFromEPSG(epsg_from_num)

                    utm_band = str((math.floor((ulx_trns + 180)/6)%60)+1)
                    if  len(utm_band) == 1:
                        utm_band = '0'+ utm_band
                    if uly_trns >= 0:
                        epsg_code = '326' + utm_band
                    else:
                        epsg_code = '327' + utm_band

                    epsg_code = int(epsg_code)
                    outSpatial_ref = osr.SpatialReference()
                    outSpatial_ref.ImportFromEPSG(epsg_code)

                    Geo_to_Proj = osr.CoordinateTransformation(inSpatial_ref, outSpatial_ref)        

                    (UlX, UlY, UlZ) = Geo_to_Proj.TransformPoint(gt[0], gt[3])
                    (LrX, LrY, LrZ) = Geo_to_Proj.TransformPoint(gt[0] + gt[1]*xRS, gt[3] + gt[5]*yRS)

                    num_of_pixels = int(round((LrX - UlX)/outPixelsize, 0))
                    outPixelsize = abs(((gt[0] + gt[1]*xRS) - ulx_trns)/num_of_pixels)
                    
                elif inCStype == 'PROJCS':
                    
                    num_of_pixels = ((gt[0] + gt[1]*xRS) - ulx_trns)/pixel_size 
                    outPixelsize = (lrx - ulx)/num_of_pixels
               
                
        if data_type == '8bit':
            #NoDataValue = 0                                                                                                                         #8BIT HERE
            proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", int((lrx-ulx)/outPixelsize), int((uly-lry)/outPixelsize), 1, gdal.GDT_Byte)        #8BIT HERE
        elif data_type == '32bit':
            proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", int((lrx-ulx)/outPixelsize), int((uly-lry)/outPixelsize), 1, gdal.GDT_Float32)     #32BIT HERE
        else:
            proj_ras = gdal.GetDriverByName( 'MEM' ).Create( "", int((lrx-ulx)/outPixelsize), int((uly-lry)/outPixelsize), 1, gdal.GDT_Float32)     #32BIT HERE

        new_gt = (ulx, outPixelsize, gt[2], uly, gt[4], -outPixelsize)
        proj_ras.SetGeoTransform(new_gt)
        proj_ras.SetProjection(to_coord.ExportToWkt())
        proj_ras_band = proj_ras.GetRasterBand(1)
        proj_ras_band.SetNoDataValue(NoDataValue)
        
        nodata_array = np.full((int((uly-lry)/outPixelsize), int((lrx-ulx)/outPixelsize)), NoDataValue)     # NoData HACK
        proj_ras_band.WriteArray(nodata_array)                                                              # NoData HACK

        if resample_alg == 'Nearest':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_NearestNeighbour) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Bilinear':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Bilinear) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Cubic':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_Cubic) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
        elif resample_alg == 'Cubic_spline':
            res = gdal.ReprojectImage(inData, proj_ras, from_coord.ExportToWkt(), to_coord.ExportToWkt(), gdal.GRA_CubicSpline) # GRA_NearestNeighbour, GRA_CubicSpline, GRA_Bilinear, GRA_Cubic, GRA_Lanczos  
            

        gdal.GetDriverByName('GTiff').CreateCopy(outRas, proj_ras)
        inData = None
        proj_ras = None
        
    #if snap != False:
    #    print('\tSnapping ...')
    #    snap_raster_prj(snap, proj_ras, outRas, to_coord.ExportToWkt())
    #    proj_ras=None
    #else:
    #    gdal.GetDriverByName('GTiff').CreateCopy(outRas, proj_ras)
    #    proj_ras=None


def get_raster_epsg(inraster):
    '''
    Function Def: get_raster_epsg(inraster)

    inraster     ~ path to a raster file [REQUIRED]
    - - - - - - - - - - - - - - - - - - - -

    returns      ~ list containing epsg code, if found, and a confidence ranking (1 to 6, with 1 being the highes confidence) .e.g ['EPSG:32734', 1]
                 ~ None, if the EPSG code could not be found.
    
    '''
    raster_epsg = wkt2epsg(get_raster_wkt(inraster), forceProj4 = False, cstypeStat = False)
    return raster_epsg

def get_vector_epsg(invector):
    '''
    Function Def: get_vector_epsg(invector)

    invector     ~ path to a vector file [REQUIRED]
    - - - - - - - - - - - - - - - - - - - -

    returns      ~ list containing epsg code, if found, and a confidence ranking (1 to 6, with 1 being the highes confidence) .e.g ['EPSG:32734', 1]
                 ~ None, if the EPSG code could not be found.
    
    '''
    vector_epsg = wkt2epsg(get_shapefile_wkt(invector), forceProj4 = False, cstypeStat = False)
    return vector_epsg



def get_pixel_size_m(inRas):

    CStype = wkt2epsg(get_raster_wkt(inRas), cstypeStat = True)   # [IN] Checks the coordinate system type
    
    
    if CStype == 'GEOGCS': # if coord is geographical
        #________________________________________________________looks for the UTM zone EPSG code (pretty schweet)

        inData = gdal.Open(inRas)
        gt = inData.GetGeoTransform()       # Geotransformation # ulx[0], xres[1], xskew[2], uly[3], yskew[4], yres[5]
        x_arr_size = inData.RasterXSize
        y_arr_size = inData.RasterYSize

        x_min = gt[0]
        x_max = gt[0] + (gt[1]*x_arr_size)

        y_min = gt[3] + gt[5]*y_ar_size
        y_max = gt[3]
        inData = None

        #x_min, x_max, y_min, y_max = vector_layer.GetExtent()   # [IN] Get the extent of the vector file (Extent coordinates will be in the format of the input vector files projection)[Layer]

        ePsG_list = get_raster_epsg(inRas)
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
        x_pixel = pixel_size                             
        y_pixel = -pixel_size

        return x_pixel, y_pixel

    else:
        inData = gdal.Open(inRas)
        gt = inData.GetGeoTransform()       # Geotransformation # ulx[0], xres[1], xskew[2], uly[3], yskew[4], yres[5]
        x_pixel = int(round(gt[1],0))
        y_pixel = int(round(gt[5],0))

        return x_pixel, y_pixel

        
        
#DEM_wgs = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_wgs.tif'
#espg, confidence = get_raster_epsg(DEM_wgs)


'''

#### TESTING THE SCRIPT VECTOR

polygon = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/square.shp'
polygonOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/squarePRJ.shp'

print('Projecting shape polygon1')
project_shp(polygon, polygonOUT, 'EPSG:102022', 'polygon')
print('Done')

singlep = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/single_poly.shp'
singlepOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/single_poly_PRJ.shp'

print('Projecting shape polygon2')
project_shp(singlep, singlepOUT, 'EPSG:102022', 'polygon')
print('Done')

point = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/TMP.shp'
pointOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/TMP_PRJ.shp'

print('Projecting shape point')
project_shp(point, pointOUT, 'EPSG:102022', 'point')
print('Done')

line = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/test_line.shp'
lineOUT = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/test_linePRJ.shp'

print('Projecting shape line')
project_shp(line, lineOUT, 'EPSG:102022', 'line')
print('Done')
#project_shp(inShape, outShape, outEpsg):







#### TESTING THE SCRIPT RASTER

gdal_out2 = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_gdal_prj_50m.tif'
gdal_out = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_gdal_prj.tif'

DEM_wgs = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_wgs.tif'
DEM_prjv1 = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_prj_AAEAC.tif'
DEM_prjv2 = '//sungis08/CGA/Clients/040_Winetech/0402017001_Climate_SDSS_Pilot/Analysis/Software_Design/Test_area/DEM_prj_UTM_Zone_34S.tif'


print('Geographic coord:')
print('EPSG: ', wkt2epsg(get_raster_wkt(DEM_wgs), cstypeStat = False))
print('______________________________________________________________________________________________________\n')

print('Projected UTM:')
print('EPSG: ', wkt2epsg(get_raster_wkt(DEM_prjv2), cstypeStat = True))
print('______________________________________________________________________________________________________\n')


print('Projected esri custom AAEAC:')
#print(get_raster_wkt(DEM_prjv1),'\n')
print('EPSG: ', wkt2epsg(get_raster_wkt(DEM_prjv1), cstypeStat = False))
print('______________________________________________________________________________________________________\n')


print('Projected gdal EPSG AAEAC:')
print('EPSG: ', wkt2epsg(get_raster_wkt(gdal_out), cstypeStat = True))
print('______________________________________________________________________________________________________\n')

      
print('Start projecting ...')
project_raster(DEM_wgs, gdal_out2, 'EPSG:102022', outPixelsize = 50)
#(inRas, outRas, outEpsg, outPixelsize=None, snap=False):
#project_raster(DEM_prjv1, gdal_out2, 'EPSG:4326')
#wkttext  = get_raster_wkt(DEM_prjv1)

print('Done')

'''

