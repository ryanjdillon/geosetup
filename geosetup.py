#/usr/bin/env python

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import datetime
import numpy.lib.recfunctions
import pyproj
from StringIO import StringIO
import geosetup.mpl_util
from geosetup.globcolour import plotglob_mapped
from geosetup.cortad import getCortad
from geosetup.griddata import invdistgis
from geosetup.griddata import invdist
from geosetup.griddata import data2raster
from geosetup.griddata.datainterp import geointerp

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
    b_dist = g.inv(west_lon,north_lat,east_lon,north_lat)[2]/2
    c_dist = g.inv(west_lon,north_lat,lon0,south_lat)[2]

    mapH = pow(pow(c_dist,2)-pow(b_dist,2),1./2)
    lat0 = g.fwd(lon0,south_lat,0,mapH/2)[1]

    # distance between max E and W longitude at most southern latitude
    mapW = g.inv(west_lon,south_lat,east_lon,south_lat)[2]

    return lon0,lat0,mapW*scale,mapH*scale

def drawsst(ax,map_object,longSST,latSST,filledSST,myalpha):
    '''
    Draw Sea surface temperature - Trond Kristiansen
    '''
    # Input arrays have to be 2D
    print "Drawing SST: max %s and min %s"%(filledSST.min(), filledSST.max())
    x2, y2 = map_object(longSST,latSST)
    levels=np.arange(2,18,0.5)

    # TODO correct issue with colormap
    if myalpha > 0.99:
        CS2 = map_object.contourf(x2,y2,filledSST,levels,
                #cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r),
                cmap=plt.cm.jet,
                extend='upper',alpha=myalpha)
    else:
        CS2 = map_object.contourf(x2,y2,filledSST,levels,
                #cmap=mpl_util.LevelColormap(levels,cmap=cm.RdYlBu_r),
                cmap=plt.cm.jet,
                extend='upper',alpha=myalpha)

        # TODO correct or remove
        #CS2 = map_object.contourf(x2,y2,filledSST,levels,
        #               cmap=mpl_util.LevelColormap(levels,cmap=cm.Greys),
        #               extend='upper',alpha=myalpha)

def find_nearest(array,value):
    '''
    TODO
    '''
    idx = (np.abs(array-value)).argmin()
    #return array[idx] # return value
    return idx # return index

def filter2bool(regexp,array):
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
    BTM_DIR = 'data/gebco/'
    GRD_SIZE = 50 #km
    GRD_LAT_START = 40.
    GRD_LAT_STOP = 50.
    GRD_LON_START = -10.
    GRD_LON_STOP = 10.

    #########################
    # Process Sighting Data #
    #########################

    data_file = sys.argv[1]

    if sys.argv[2]:
        rows_to_skip = int(sys.argv[2])
    else:
        rows_to_skip = 0

    # TODO double check if following is necessary
    with open(data_file) as fh:
        io = StringIO(fh.read().replace(',', '\t'))

    # Define names and data types for sighting data
    record_types = np.dtype([
                    ('vessel',str,1),           #00 - Vessel ID
                    ('dates',str,6),            #01 - Date
                    ('times',str,6),            #02 - Time (local?)
                    ('lat',float),              #03 - latitude dec
                    ('lon',float),              #04 - longitude dec -1
                    ('beafort',str,2),          #05 - beafort scale
                    ('weather',int),            #06 - weather code
                    ('visibility',int),         #07 - visibility code
                    ('effort_sec',float),       #08 - seconds on effort
                    ('effort_nmil',float),      #09 - n miles on effort
                    ('lat0',float),             #10 - lat of sight start
                    ('lon0',float),             #11 - lon of sight start
                    ('num_observers',int),      #12 - number of observers
                    ('species',str,6),          #13 - species codes
                    ('num_animals',int),        #14 - number observed
                    ('sighting',int),           #15 - Boolean sight code
                    ('rdist',float),            #16 - distance to sight
                    ('angle',float),            #17 - angel from ship
                    ('block',str,2),       #18 - cruise block
                    ('leg',int),                #19 - cruise leg code
                    ('observations',str,25),    #20 - cruise obs codes
                    ])

    # Import data to structured array
    data = np.genfromtxt(io,dtype=record_types,delimiter='\t',skip_header=rows_to_skip)

    # Correct longitude values
    data['lon'] = data['lon']*(-1)

    # Create array of unique survey block IDs
    blocks = np.unique(data['block'])

    # Get date information for subsampling environment data
    dates = list()
    for data_date,data_time in zip(data['dates'],data['times']):
        data_datetime = data_date+data_time
        dates.append(datetime.datetime.strptime(data_datetime, '%y%m%d%H%M%S'))

    data_start = min(dates)
    data_end = max(dates)

    # Print a summary of geo data 
    print '\nSighting Data Information'
    print '-------------------------------------------'
    print 'Data Path: '+data_file
    print 'First sighting: ',data_start
    print 'Last sighting: ',data_end


    #######################################
    # Calculate Sightings per unit effort #
    #######################################

    # 'BM' Blue Whale (Balaenoptera musculus)
    # 'BP' Fin Whale (Balaenoptera physalus)
    # 'BB' Sei Whale (Balaenoptera borealis)
    # 'BA' Minke whales (Balaenoptera acutorostrata)
    # 'MN' Humpback Whale (Megaptera novaeangliae)

    # TODO calculate distance between start points and sighting point and compare
#    g = pyproj.Geod(ellps='WGS84') # Use WGS84 ellipsoid
#    f_azimuth, b_azimuth, dist = g.inv(data['lon'],data['lat'],data['lon0'],data['lat0'])

    # get create array of indexes for minke whales
    # TODO verify the effort calc is correct
    minke_idx = np.where(filter2bool('BA',data['species'])==True)[0]
    minke_effort = data['effort_nmil'][minke_idx]
    minke_lat = data['lat'][minke_idx]
    minke_lon = data['lon'][minke_idx]

    minke_proj = data2raster.GeoPoint(minke_lon,minke_lat,minke_effort)
    minke_x, minke_y = minke_proj.transform_point()

    ZI = invdist.invDist(minke_lat,minke_lon, minke_effort)

    XI, YI = np.meshgrid(minke_x, minke_y)
    n = plt.normalize(0.0, 1000.0)
    plt.subplot(1, 1, 1)
    plt.pcolor(XI, YI, ZI)
#    plt.scatter(xv, yv, 100, values)
    plt.colorbar()
    plt.show()

    # Create Effort Gtiff
    #spueGeopoint = data2raster.GeoPoint(data['lon'],data['lat'],data['spue'])
    #spueGeopoint.create_raster(filename="spue.tiff",output_format="GTiff")

    #TODO remove following
#    last_idx = 0
#    idx_pos = 0
#    minke_effort = np.zeros_like(minke_idx, dtype=float)
#    for idx in minke_idx:
#        minke_effort[idx_pos] = dist[last_idx:(idx+1)].sum()
#        last_idx = idx
#        idx_pos = idx_pos+1

    # generate list of indexes where effort was greater than zero
#    effort_idx = np.where(data['effort_sec']*data['effort_nmil'] != 0)

#    spue = data['num_animals'][effort_idx]/(data['effort_sec'][effort_idx]*data['effort_nmil'][effort_idx])

    # append spue calculations to structured array dataset
#    data = numpy.lib.recfunctions.append_fields(data,'spue',data=spue)

    #######################
    # Process CorTAD Data #
    #######################

    ref_date=datetime.datetime(1980,12,31,12,0,0)
    days = 60.*60.*24. # sec*min*hr

    # calculate days from ref date to first sighting
    time_start = int(round((data_start - ref_date).total_seconds()/days))
    # calculate days from ref date to last sighting
    time_end = int(round((data_end - ref_date).total_seconds()/days))

    filledSST = getCortad.extractCORTADSST("North Sea",time_start,time_end)

    lonSST2D, latSST2D, lonSST, latSST = getCortad.extractCoRTADLongLat()

    # Print sighting data informtion
    print 'Sighting Period: ', time_start, time_end, data_end-data_start

    # Create SST Gtiff
    #sstGeopoint = data2raster.GeoPoint(data['lon'],data['lat'],data['spue'])
    #sstGeopoint.create_raster(filename="spue.tiff",output_format="GTiff")

    #################################
    # Process Globcolour Chl-a Data #
    #################################

    # Extract Chl-a data
    chla_lons, chla_lats, chla_vals = plotglob_mapped.getMappedGlobcolour(PROJ_DIR+CHL_DIR,0.,50.,30.,60.,'2007-08-01','2007-08-30')
    chla_vals = np.ravel(chla_vals)

    # Create Chla Gtiff
    chlaGeopoint = data2raster.GeoPoint(chla_lons,chla_lats,chla_vals)
    chlaGeopoint.create_raster(filename="chla.tiff",output_format="GTiff",
cell_width_meters = 50000, cell_height_meters = 50000)

    ############################
    # Process Bathymetric Data #
    ############################

    # Create Bathy Gtiff
    #bathyGeopoint = data2raster.GeoPoint(data['lon'],data['lat'],data['spue'])
    #bathyGeopoint.create_raster(filename="bathy.tiff",output_format="GTiff")

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
    drawsst(ax,m,lonSST2D,latSST2D,filledSST,0.25)

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
