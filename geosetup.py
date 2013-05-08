#/usr/bin/env python

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import sys, os, errno
import re
import math
import datetime
import numpy.lib.recfunctions
import pyproj
from StringIO import StringIO
import geosetup.mpl_util
#TODO make module naming consistent
from geosetup.globcolour import globcolour
from geosetup.pathfinder import pathfinder
from geosetup.gebco import gebco
from geosetup.interpolate import invdistgis, invdist
from geosetup.interpolate.datainterp import geointerp
from geosetup.writefile import data2raster
from geosetup.sightsurvey import sightsurvey

# prevent creation of .pyc files
sys.dont_write_bytecode = True

#############
# Functions #
#############

def plotSizedData(map_obj,lons,lats,values,symbol,min_size,max_size,lines=False):
    '''
    Plot data with varying sizes
    '''
    proj_x,proj_y = m(lons,lats)

    # Calculate marker sizes using y=mx+b
    # where y = marker size and x = data value
    slope = (max_size-min_size)/(max(values)-min(values))
    intercept = min_size-(slope*min(values))

    for x, y, val in zip(proj_x, proj_y, values):
        msize = (slope*val)+intercept
        map_obj.plot(x, y, symbol, markersize=msize)

def plotLines(map_obj,lons,lats,lw=1.0,color='k'):
    '''
    Draw a line between each set of coordinates
    '''
    for i in range(len(lons)-1):
        map_obj.drawgreatcircle(lons[i],lats[i],lons[i+1],lats[i+1],linewidth=1.0,color='k')

def centerMap(lons,lats,scale):
    '''
    Set range of map. Assumes -90 < Lat < 90 and -180 < Lon < 180, and
    latitude and logitude are in decimal degrees
    '''
    north_lat = max(lats)
    south_lat = min(lats)
    west_lon = max(lons)
    east_lon = min(lons)

    # find center of data
    # average between max and min longitude
    lon0 = ((west_lon-east_lon)/2.0)+east_lon

    # define ellipsoid object for distance measurements
    g = pyproj.Geod(ellps='WGS84') # Use WGS84 ellipsoid TODO make variable
    earth_radius = g.a # earth's radius in meters

    # Use pythagorean theorom to determine height of plot
    # divide b_dist by 2 to get width of triangle from center to edge of data area
    # inv returns [0]forward azimuth, [1]back azimuth, [2]distance between

    # a_dist = the height of the map (i.e. mapH)
    b_dist = g.inv(west_lon, north_lat, east_lon, north_lat)[2]/2
    c_dist = g.inv(west_lon, north_lat, lon0, south_lat)[2]

    mapH = pow(pow(c_dist,2)-pow(b_dist,2),1./2)
    lat0 = g.fwd(lon0,south_lat,0,mapH/2)[1]

    # distance between max E and W longitude at most southern latitude
    mapW = g.inv(west_lon, south_lat, east_lon, south_lat)[2]

    return lon0, lat0, mapW*scale, mapH*scale

def drawsst(ax, map_object, longSST, latSST, filledSST, myalpha):
    '''
    Draw Sea surface temperature - Trond Kristiansen
    '''
    # Input arrays have to be 2D
    print "Drawing SST: max %s and min %s"%(filledSST.min(), filledSST.max())
    x2, y2 = map_object(longSST,latSST)
    levels=np.arange(2,18,0.5)

    # TODO correct issue with colormap
    if myalpha > 0.99:
        CS2 = map_object.contourf(x2, y2, filledSST, levels,
                #cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r),
                cmap = plt.cm.jet,
                extend = 'upper', alpha=myalpha)
    else:
        CS2 = map_object.contourf(x2, y2, filledSST, levels,
                #cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r),
                cmap = plt.cm.jet,
                extend = 'upper', alpha=myalpha)

        # TODO correct or remove
        #CS2 = map_object.contourf(x2,y2,filledSST,levels,
        #               cmap=mpl_util.LevelColormap(levels,cmap=cm.Greys),
        #               extend='upper',alpha=myalpha)

def find_nearest(array, value):
    '''
    TODO
    '''
    idx = (np.abs(array - value)).argmin()
    #return array[idx] # return value
    return idx # return index

def filter2bool(regexp, array):
    '''
    Create array of positive boolean where elements match regex 
    '''
    return np.array([bool(re.search(regexp, element)) for element in array])

if __name__ == '__main__':

    #####################
    # Commandline Usage #
    #####################

    if len(sys.argv) < 3:
        print >>sys.stderr,'\nUsage:',sys.argv[0],'<datafile> <#rows to skip>\n'
        sys.exit(1)

    ############################
    # Configuration parameters #
    ############################

    PROJ_DIR = '/home/ryan/Desktop/asf-fellowship/code/geosetup/'
    SST_DIR = 'data/pathfinder/'
    CHL_DIR = 'data/globcolour/'
    BTM_DIR = 'data/gebco/gridone.nc'
    SIGHT_DATA = sys.argv[1] # 'data/survey/na07.tab'
    OUT_DIR = 'output/'
    #TODO incorporate regrid size
    GRD_SIZE = 50 #km or deg?

    # Test that output directory exists, if not create it
    try:
        os.makedirs(OUT_DIR)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

    #########################
    # Process Sighting Data #
    #########################

    # Get sight surveying data and geographical and time bounds for data
    data = sightsurvey.getData(SIGHT_DATA, skip_rows=0)

    # Get date information for subsampling environment data
    datetimes = list()
    for data_date,data_time in zip(data['dates'],data['times']):
        data_datetime = data_date+data_time
        datetimes.append(datetime.datetime.strptime(data_datetime,'%y%m%d%H%M%S'))

    LON_START = min(data['lon'])
    LON_END = max(data['lon'])
    LAT_START = min(data['lat'])
    LAT_END = max(data['lat'])
    TIME_START = min(datetimes)
    TIME_END = max(datetimes)

    # Print sighting data informtion
    print 'Sighting Period: ', TIME_START, TIME_END, TIME_END - TIME_START
    ###############
    # Create Grid #
    ###############

    grid_lat_start =  math.floor(LAT_START)
    grid_lon_start =  math.floor(LON_START)
    grid_lat_end =  math.ceil(LAT_END)
    grid_lon_end =  math.ceil(LON_END)
    grid_lats = np.linspace(GRD_SIZE, 180, 180/GRD_SIZE)
    grid_lons = np.linspace(GRD_SIZE, 180, 180/GRD_SIZE)

    #######################################
    # Calculate Sightings per unit effort #
    #######################################
    # 'BM' Blue Whale (Balaenoptera musculus)
    # 'BP' Fin Whale (Balaenoptera physalus)
    # 'BB' Sei Whale (Balaenoptera borealis)
    # 'BA' Minke whales (Balaenoptera acutorostrata)
    # 'MN' Humpback Whale (Megaptera novaeangliae)

    # Create array of indexes for minke whales
    # TODO verify the effort calc is correct
    minke_idx = np.where(filter2bool('BA',data['species'])==True)[0]
    minke_effort = data['effort_nmil'][minke_idx]
    minke_lat = data['lat'][minke_idx]
    minke_lon = data['lon'][minke_idx]

    # Create Effort Gtiff
    effortGeopoint = data2raster.GeoPoint(minke_lon, minke_lat, minke_effort)
    effortGeopoint.create_raster(filename = OUT_DIR + "effort.tiff", output_format="GTiff")

    # Interpolate / Plot Minke effort
    minke_x, minke_y = effortGeopoint.transform_point()
    ZI = invdist.invDist(minke_lat,minke_lon, minke_effort)

    XI, YI = np.meshgrid(minke_x, minke_y)
    n = plt.normalize(0.0, 1000.0)
    plt.subplot(1, 1, 1)
    plt.pcolor(XI, YI, ZI)
    #plt.scatter(xv, yv, 100, values)
    plt.colorbar()
    plt.show()

    #######################
    # Process CorTAD Data # TODO revue / remove
    #######################

#    # Get cortad SST within date period
#    # TODO modify method to subset lat/lon
#    filledSST = cortad.extractCORTADSST("North Sea", TIME_START, TIME_END)
#    lonSST2D, latSST2D, sst_lon, sst_lat = cortad.extractCoRTADLongLat()
#    sst_lon, sst_lat = np.meshgrid(sst_lon,sst_lat)
#    sst_lon = np.ravel(sst_lon)
#    sst_lat = np.ravel(sst_lat)
#    filledSST_flat = np.ravel(filledSST)

#    # Create SST Gtiff
#    sstGeopoint = data2raster.GeoPoint(sst_lon, sst_lat, filledSST_flat)
#    sstGeopoint.create_raster(filename = OUT_DIR + "sst.tiff", output_format="GTiff")

    ###############################
    # Process Pathfinder SST Data #
    ###############################

    # Extract Chl-a datai
    # getnetcdfdata(data_dir, nc_var_name, min_lon, max_lon, min_lat, max_lat, data_time_start, data    _time_end):
    sst_lons, sst_lats, sst_vals = pathfinder.getnetcdfdata(PROJ_DIR + SST_DIR, 
                                                              'sea_surface_temperature',
                                                               LON_START, LON_END,
                                                               LAT_START, LAT_END,
                                                               TIME_START, TIME_END)
    # Create SST Gtiff
    sstGeopoint = data2raster.GeoPoint(sst_lons, sst_lats, sst_vals)
    sstGeopoint.create_raster(filename = OUT_DIR + "sst.tiff", output_format="GTiff",
                               cell_width_meters = 5000, cell_height_meters = 5000)

    #################################
    # Process Globcolour Chl-a Data #
    #################################

    # Extract Chl-a data
    chla_lons, chla_lats, chla_vals = globcolour.getMappedGlob(PROJ_DIR+CHL_DIR,
                                                               LON_START, LON_END,
                                                               LAT_START, LAT_END,
                                                               TIME_START, TIME_END)
    chla_vals = np.ravel(chla_vals)

    # Create Chla Gtiff
    chlaGeopoint = data2raster.GeoPoint(chla_lons,chla_lats,chla_vals)
    chlaGeopoint.create_raster(filename = OUT_DIR + "chla.tiff", output_format="GTiff",
                               cell_width_meters = 5000, cell_height_meters = 5000)

    ############################
    # Process Bathymetric Data #
    ############################

    # Extract bathymetric data
    bathy_lons, bathy_lats, bathy_z = gebco.getGebcoData(BTM_DIR, LON_START, LON_END,
                                                                  LAT_START, LAT_END)
    # Create Bathy Gtiff
    bathyGeopoint = data2raster.GeoPoint(bathy_lons,bathy_lats,bathy_z)
    bathyGeopoint.create_raster(filename = OUT_DIR + "bathy.tiff",output_format="GTiff")

    ###############
    # Create Plot #
    ###############

    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111)

    # Calculate map's center lat and lon from sampling data
    lon0center,lat0center,mapWidth,mapHeight = centerMap(data['lon'],data['lat'],1.1)

    # Lambert Conformal Projection Plot
    # lat_1 is first standard parallel.
    # lat_2 is second standard parallel (defaults to lat_1).
    # lon_0,lat_0 is central point.
    # TODO check that following isn't more accurate from pyProj
    # rsphere=(6378137.00,6356752.3142) specifies WGS4 ellipsoid
    # area_thresh=1000 means don't plot coastline features less
    # than 1000 km^2 in area.
    m = Basemap(width=mapWidth,height=mapHeight,
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',area_thresh=1000.,projection='lcc',\
                lat_1=45.,lat_2=55.,\
                lat_0=lat0center,lon_0=lon0center)

    # Plot STT
    #drawsst(ax,m,lonSST2D,latSST2D,filledSST,0.25)

    # Plot Chl-a
    #TODO
    #drawchla()

    # Plot Depth
    #TODO
    #drawbottom()

    # Draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,0,0,1], fontsize=10)
    m.drawmapboundary(fill_color='aqua')
    m.drawcoastlines(linewidth=0.2)
    m.fillcontinents(color='white', lake_color='aqua')

    # Plot data points
    plotSizedData(m,minke_lon,minke_lat,minke_effort,'ro',4,20)
    plotLines(m,data['lon'],data['lat'],lw=1.0,color='k')
    #x, y = m(data['lon'],data['lat'])
    #m.scatter(x,y,2,marker='o',color='k')

    # Print Plot Information
    print '\nPlot Information'
    print '-------------------------------------------'
    print "Map lat/lon center: ",lat0center,lon0center
    print "Map height/width: ",mapWidth,mapHeight
    # TODO write metadata to image file
    # http://stackoverflow.com/questions/10532614/can-matplotlib-add-metadata-to-saved-figures
    #plt.title("Example")
    #plotfile='SST_northsea_'+str(currentDate)+'.png'
    #print "Saving map to file %s"%(plotfile)
    #plt.savefig(plotfile)

    plt.show()
