from osgeo import gdal, osr
import numpy
import math

#############
# Functions #
#############

def creategrid(min_lon, max_lon, min_lat, max_lat, cell_size_deg, mesh=False):
    '''Output grid within geobounds and specifice cell size

       cell_size_deg should be in decimal degrees'''

    min_lon = math.floor(min_lon)
    max_lon = math.ceil(max_lon)
    min_lat = math.floor(min_lat)
    max_lat = math.ceil(max_lat)

    lon_num = (max_lon - min_lon)/cell_size_deg
    lat_num = (max_lat - min_lat)/cell_size_deg

    grid_lons = np.zeros(lon_num) # fill with lon_min
    grid_lats = np.zeros(lat_num) # fill with lon_max
    grid_lons = grid_lons + (np.assary(range(lon_num))*cell_size_deg)
    grid_lats = grid_lats + (np.assary(range(lat_num))*cell_size_deg)

    grid_lons, grid_lats = np.meshgrid(grid_lons, grid_lats)
    grid_lons = np.ravel(grid_lons)
    grid_lats = np.ravel(grid_lats)

    #if mesh = True:
    #    grid_lons = grid_lons
    #    grid_lats = grid_lats

    return grid_lons, grid_lats

geo = [-180, 1./6, 0, 90., 0., -1./6]
nx = 360*6
ny = 180*6
driver = gdal.GetDriverByName ("GTiff")
dst_ds = driver.Create ( "/tmp/world_grid.tif",nx, ny, 1, gdal.GDT_Int32)
dst_ds.SetGeoTransform ( geo )
proj = osr.SpatialReference()
proj.ImportFromEPSG ( 4326)
dst_ds.SetProjection ( proj.ExportToWkt() )
dst_ds.GetRasterBand(1).WriteArray ( numpy.zeros ( (ny, nx) ))
dst_ds = None
