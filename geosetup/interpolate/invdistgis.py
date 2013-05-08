from math import pow
from math import sqrt
import numpy as np
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from osgeo.gdalconst import *
import sys
''' http://geoexamples.blogspot.com/2012/05/creating-grid-from-scattered-data-using '''

def pointValue(x,y,power,smoothing,xv,yv,values):
    ''' interpolation algorithm re-written from gdal libraries'''
    nominator=0
    denominator=0
    for i in range(0,len(values)):
        dist = sqrt((x-xv[i])*(x-xv[i])+(y-yv[i])*(y-yv[i])+smoothing*smoothing);
        #If the point is really close to one of the data points, return the data point value to avoid singularities
        if(dist<0.0000000001):
            return values[i]
        nominator=nominator+(values[i]/pow(dist,power))
        denominator=denominator+(1/pow(dist,power))
    #Return NODATA if the denominator is zero
    if denominator > 0:
        value = nominator/denominator
    else:
        value = -9999
    return value

def invDist(xv,yv,values,geotransform,proj,xSize,ySize,power,smoothing,driverName,outFile):
    #Transform geographic coordinates to pixels
    for i in range(0,len(xv)):
        xv[i] = (xv[i]-geotransform[0])/geotransform[1]
    for i in range(0,len(yv)):
        yv[i] = (yv[i]-geotransform[3])/geotransform[5]
    #Creating the file
    driver = gdal.GetDriverByName( driverName )
#    ds = driver.Create( outFile, xSize, ySize, 1, gdal.GDT_Float32)
    ds = driver.Create( outFile, int(xSize), int(ySize), int(1.0), gdal.GDT_Float32)
    if proj is not None:
        ds.SetProjection(proj.ExportToWkt())
    ds.SetGeoTransform(geotransform)
    valuesGrid = np.zeros((ySize,xSize))
    #Getting the interpolated values
    for x in range(0,xSize):
        for y in range(0,ySize):
            valuesGrid[y][x] = pointValue(x,y,power,smoothing,xv,yv,values)

    ds.GetRasterBand(1).WriteArray(valuesGrid)
    ds = None
    return valuesGrid

def readPoints(dataFile,Zfield='Z'):
    data = {}
    xv=[]
    yv=[]
    values=[]
    ds = ogr.Open(dataFile)
    if ds is None:
        raise Exception('Could not open ' + dataFile)


    layer = ds.GetLayer()
    proj = layer.GetSpatialRef()
    extent = layer.GetExtent()

    feature = layer.GetNextFeature()
    if feature.GetFieldIndex(zField) == -1:
        raise Exception('zField is not valid: ' + zField)

    while feature:
        geometry = feature.GetGeometryRef()
        xv.append(geometry.GetX())
        yv.append(geometry.GetY())
        values.append(feature.GetField(zField))

        feature = layer.GetNextFeature()
    data['extent'] = extent
    data['xv']=xv
    data['yv']=yv
    data['values']=values
    data['proj'] = proj
    ds = None
    return data

def usage():
    print "Usage: python indvdistgis.py [-zfield field_name] [-a_srs srs_def] [-out_format out_format] [-outsize xsize ysize] [-txe xmin xmax] [-tye ymin ymax] data_file out_file"
    sys.exit(1)

if __name__ == "__main__":
    #Setting the default values
    power=2
    smoothing=0

    zField='Z'
    dataFile=None
    outFile=None
    driverName='EHdr' #'GTiff'
    proj=None

    geotransform = None
    xMin=None
    xMax=None
    yMin=None
    yMax=None
    xSize=100
    ySize=100

    #Parsing the command line
    argv = gdal.GeneralCmdLineProcessor( sys.argv )

    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg == '-out_format':
            driverName = argv[i+1]
            driverName = driverName.replace("'","")
            driverName = driverName.replace('"','')
            i = i + 1
        elif arg == '-zfield':
            zField = argv[i+1]
            zField = zField.replace("'","")
            zField = zField.replace('"','')
            i = i + 1
        elif arg == '-a_srs':
            proj = argv[i+1]
            proj = proj.replace("'","")
            proj = proj.replace('"','')
            i = i + 1
        elif arg == '-outsize':
            xSize = argv[i+1]
            xSize = xSize.replace("'","")
            xSize = int(xSize.replace('"',''))
            ySize = argv[i+2]
            ySize = ySize.replace("'","")
            ySize = int(ySize.replace('"',''))
            i = i + 2
        elif arg == '-txe':
            xMin = argv[i+1]
            xMin = xMin.replace("'","")
            xMin = float(xMin.replace('"',''))
            xMax = argv[i+2]
            xMax = xMax.replace("'","")
            xMax = float(xMax.replace('"',''))
            i = i + 2
        elif arg == '-tye':
            yMin = argv[i+1]
            yMin = yMin.replace("'","")
            yMin = float(yMin.replace('"',''))
            yMax = argv[i+2]
            yMax = yMax.replace("'","")
            yMax = float(yMax.replace('"',''))
            i = i + 2
        elif dataFile is None:
            dataFile = arg
            dataFile = dataFile.replace("'","")
            dataFile = dataFile.replace('"','')
        elif outFile is None:
            outFile = arg
            outFile = outFile.replace("'","")
            outFile = outFile.replace('"','')
        i = i + 1

    if dataFile is None or outFile is None:
        usage()
    try:
        data = readPoints(dataFile,zField)
    except Exception,ex:
        print ex
        sys.exit(0)
    if xMin is None:
        xMin=data['extent'][0]
        xMax=data['extent'][1]
    if yMin is None:
        yMin=data['extent'][2]
        yMax=data['extent'][3]

    geotransform=[]
    geotransform.append(xMin)
    geotransform.append((xMax-xMin)/xSize)
    geotransform.append(0)
    geotransform.append(yMax)
    geotransform.append(0)
    geotransform.append((yMin-yMax)/ySize)

    if proj is None:
        proj = data['proj']
        print proj
        print data['extent']
        print data['xv']
        print data['yv']
    else:
        try:
            proj = osr.SpatialReference()
            proj.SetFromUserInput(str(proj))
        except Exception,ex:
            print ex
            sys.exit(0)
    #Creating the interpolation function and populating the output matrix value
    try:
        ZI = invDist(data['xv'],data['yv'],data['values'],geotransform,proj,xSize,ySize,power,smoothing,driverName,outFile)
    except Exception,ex:
        print ex
        sys.exit(0)
