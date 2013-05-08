import numpy as np
import gdal
import os
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
from osgeo.gdalconst import *
import pyproj
import scipy.sparse
import scipy
gdal.AllRegister()
gdal.UseExceptions()

'''http://monkut.webfactional.com/blog/archive/2012/5/2/understanding-raster-basic-gis-concepts-and-the-python-gdal-library/'''

#############
# Functions #
#############

def get_iterable_extent(*args):
    '''Returns list of minimum and maximum from lists/array input'''
    iterable_extent = list()
    for iter_object in args:
        iterable_extent.append(min(iter_object))
        iterable_extent.append(max(iter_object))
    return iterable_extent

class GeoPoint:
    '''Tranform coordinate system, projection, and write geospatial data'''

    """lon/lat values in WGS84"""
    # LatLon with WGS84 datum used by GPS units and Google Earth
    wgs84 = pyproj.Proj("+init=EPSG:4326")
    # Lambert Conformal Conical (LCC)
    lcc = pyproj.Proj("+init=EPSG:3034")
    # WGS84 Web Mercator (Auxillary Sphere; aka EPSG:900913)
    web_mercator = pyproj.Proj("+init=EPSG:3857")
    # TODO add useful projections, or hopefully make automatic

    def __init__(self, x, y, vals,inproj=wgs84, outproj=web_mercator,
cell_width_meters=50.,cell_height_meters=50.):
        self.x = x
        self.y = y
        self.vals = vals
        self.inproj = inproj
        self.outproj = outproj
        self.cellw = cell_width_meters
        self.cellh = cell_height_meters

    def __del__(self):
        class_name = self.__class__.__name__
        print class_name, "destroyed"

    # TODO update so following two funtions take either point or list of points
    def transform_point(self, x=None, y=None):
        if x is None:
            x = self.x
        if y is None:
            y = self.y
        return pyproj.transform(self.inproj,self.outproj,x,y)


    def get_raster_size(self, point_x, point_y, cell_width_meters, cell_height_meters):
        """Determine the number of rows/columns given the bounds of the point
        data and the desired cell size"""
        # TODO convert points to 2D projection first

        cols = int(((max(point_x) - min(point_x)) / cell_width_meters)+1)
        rows = int(((max(point_y) - min(point_y)) / abs(cell_height_meters))+1)

        print 'cols: ',cols
        print 'rows: ',rows

        return cols, rows


    def create_geotransform(self,x_rotation=0,y_rotation=0):
        '''Create geotransformation for converting 2D projected point from
        pixels and inverse geotransformation to pixels (what I want, but need the
        geotransform first).'''

        # TODO make sure this works with negative & positive lons
        top_left_x, top_left_y = self.transform_point(min(self.x),
                                                      max(self.y))

        lower_right_x, lower_right_y = self.transform_point(max(self.x),
                                                            min(self.y))
        # GeoTransform parameters
        # --> need to know the area that will be covered to define the geo tranform
        # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        # NOTE: cell height must be negative (-) to apply image space to map
        geotransform = [top_left_x, self.cellh, x_rotation,
                         top_left_y, y_rotation, -self.cellh]

        # for mapping lat/lon to pixel
        success, inverse_geotransform = gdal.InvGeoTransform(geotransform)
        if not success:
            print 'gdal.InvGeoTransform(geotransform) failed!'
            sys.exit(1)

        return geotransform, inverse_geotransform


    def point_to_pixel(self, point_x, point_y, inverse_geotransform):
        """Translates points from input projection (using the inverse
        transformation of the output projection) to the grid pixel coordinates in data
        array (zero start)"""

        # apply inverse geotranformation to convert to pixels
        pixel_x, pixel_y = gdal.ApplyGeoTransform(inverse_geotransform,
                                                  point_x, point_y)

        return pixel_x, pixel_y

    def pixel_to_point(self, pixel_x, pixel_y, geotransform):
        """Translates grid pixels coordinates to output projection points 
        (using the geotransformation of the output projection)"""

        point_x, point_y = gdal.ApplyGeoTransform(geotransform, pixel_x, pixel_y)
        return point_x, point_y

    def create_raster(self, in_x=None, in_y=None, filename="data2raster.tiff",
output_format="GTiff", cell_width_meters=1000., cell_height_meters=1000.):
        '''Create raster image of data using gdal bindings'''
        # if coords not provided, use default values from object
        if in_x is None:
            in_x = self.x
        if in_y is None:
            in_y = self.y

        # create empty raster
        current_dir = os.getcwd()+'/'
        driver = gdal.GetDriverByName(output_format)
        number_of_bands = 1
        band_type = gdal.GDT_Float32
        NULL_VALUE = 0
        self.cellw = cell_width_meters
        self.cellh = cell_height_meters

        geotransform, inverse_geotransform = self.create_geotransform()

        # convert points to projected format for inverse geotransform
        # conversion to pixels
        points_x, points_y = self.transform_point(in_x,in_y)
        cols, rows = self.get_raster_size(points_x, points_y, cell_width_meters, cell_height_meters)
        pixel_x = list()
        pixel_y = list()
        for point_x, point_y in zip(points_x,points_y):
            # apply value to array
            x, y = self.point_to_pixel(point_x, point_y, inverse_geotransform)
            pixel_x.append(x)
            pixel_y.append(y)

        dataset = driver.Create(current_dir+filename, cols, rows, number_of_bands, band_type)

        # Set geographic coordinate system to handle lat/lon
        srs = osr.SpatialReference()
        srs.SetWellKnownGeogCS("WGS84") #TODO make this a variable/generalized
        dataset.SetGeoTransform(geotransform)
        dataset.SetProjection(srs.ExportToWkt())

        #NOTE Reverse order of point components for the sparse matrix
        data = scipy.sparse.csr_matrix((self.vals,(pixel_y,pixel_x)),dtype=float)

        # get the empty raster data array
        band = dataset.GetRasterBand(1) # 1 == band index value

        # iterate data writing for each row in data sparse array
        offset = 1 # i.e. the number of rows of data array to write with each iteration
        for i in range(data.shape[0]):
            data_row = data[i,:].todense() # output row of sparse array as standard array
            band.WriteArray(data_row,0,offset*i)
        band.SetNoDataValue(NULL_VALUE)
        band.FlushCache()

        # set dataset to None to "close" file
        dataset = None

        #END FUNCTION - NO RETURN VALUE


#################
# Main Function #
#################

if __name__ == '__main__':

    # example coordinates, with function test
    lat = [45.3,50.2,47.4,80.1]
    lon = [134.6,136.2,136.9,0.5]
    val = [3,6,2,8]

    #generate some random lats and lons
#    import random
#    lats = [random.uniform(45,75) for r in xrange(500)]
#    lons = [random.uniform(-2,65) for r in xrange(500)]
#    vals = [random.uniform(1,15) for r in xrange(500)]

    lats = np.random.uniform(45,75,500)
    lons = np.random.uniform(-2,65,500)
    vals = np.random.uniform(1,25,500)
    print 'lats ',lats.shape
    #test data interpolation
    import datainterp
    lats, lons, vals = datainterp.geointerp(lats,lons,vals,2,mesh=False)
    print 'lons ',lons.shape
    geo_obj =  GeoPoint(x=lons,y=lats,vals=vals)
    geo_obj.create_raster()
