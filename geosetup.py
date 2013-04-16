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
from geosetup.cortad import getCortad
from geosetup.griddata import invdistgis
from geosetup.griddata import data2raster
from geosetup.griddata.datainterp import geointerp

# prevent creation of .pyc files
sys.dont_write_bytecode = True

#############
# Functions #
#############

def plotSizedData(lats,lons,values,plot_symbol,min_marker_size,max_marker_size):
    ''' Plot data with varying sizes '''

    # Calculate marker sizes using y=mx+b
    # where y = marker size and x = data value
    slope = (max_marker_size-min_marker_size)/(max(values)-min(values))
    intercept = min_marker_size-(slope*min(values))

    for lon, lat, val in zip(lons, lats, values):
        msize = (slope*val)+intercept
        x,y = m(lon,lat)
        m.plot(x, y, plot_symbol, markersize=msize)

def centerMap(lats,lons,scale):
    ''' Set range of map. Assumes -90 < Lat < 90 and -180 < Lon < 180, and
    latitude and logitude are in decimal degrees'''

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

    return lat0,lon0,mapW*scale,mapH*scale

def drawSST(ax,map_object,longSST,latSST,filledSST,myalpha):
    ''' Draw Sea surface temperature - Trond Kristiansen '''

    #fig = plt.figure(figsize=(12,12))
    # Input arrays has to be 2D
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

        #CS2 = map_object.contourf(x2,y2,filledSST,levels,
        #               cmap=mpl_util.LevelColormap(levels,cmap=cm.Greys),
        #               extend='upper',alpha=myalpha)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    #return array[idx] # return value
    return idx # return index

if __name__ == '__main__':
    #####################
    # Commandline Usage #
    #####################

    if len(sys.argv) < 3:
        print >>sys.stderr,'\nUsage:',sys.argv[0],'<datafile> <#rows to skip>\n'
        sys.exit(1)

    ##############
    # Setup Data #
    ##############

    data_file = sys.argv[1]
    if sys.argv[2]:
        rows_to_skip = int(sys.argv[2])
    else:
        rows_to_skip = 0

    with open(data_file) as fh:
        io = StringIO(fh.read().replace(',', '\t'))

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
                    ('blocktrack',str,2),       #18 - cruise block
                    ('leg',int),                #19 - cruise leg code
                    ('observations',str,25),    #20 - cruise obs codes
                    ])

    data = np.genfromtxt(io,dtype=record_types,delimiter='\t',skip_header=rows_to_skip)

    #######################################
    # Calculate Sightings per unit effort #
    #######################################

    # 'BM' Blue Whale (Balaenoptera musculus)
    # 'BP' Fin Whale (Balaenoptera physalus)
    # 'BB' Sei Whale (Balaenoptera borealis)
    # 'BA' Minke whales (Balaenoptera acutorostrata)
    # 'MN' Humpback Whale (Megaptera novaeangliae)

    def filter2bool(regexp,array):
        '''
        Create array of positive boolean where elements match regex 
        '''
        return np.array([bool(re.search(regexp, element)) for element in array])

    # TODO calculate distance between start points and sighting point and compare
    g = pyproj.Geod(ellps='WGS84') # Use WGS84 ellipsoid
    f_azimuth, b_azimuth, dist = g.inv(data['lon'],data['lat'],data['lon0'],data['lat0'])
    # get create array of indexs for minke whales

    #TODO verify the effort calc is correct
    minke_idx = np.where(filter2bool('BA',data['species'])==True)[0]

    last_idx = 0
    idx_pos = 0
    minke_effort = np.zeros_like(minke_idx, dtype=float)
    for idx in minke_idx:
        minke_effort[idx_pos] = dist[last_idx:(idx+1)].sum()
        last_idx = idx
        idx_pos = idx_pos+1

    # generate list of indexes where effort was greater than zero
    effort_idx = np.where(data['effort_sec']*data['effort_nmil'] != 0)

    spue = data['num_animals'][effort_idx]/(data['effort_sec'][effort_idx]*data['effort_nmil'][effort_idx])

    # append spue calculations to structured array dataset
    data = numpy.lib.recfunctions.append_fields(data,'spue',data=spue)

    # Correct longitude direction
    data['lon'] = data['lon']*(-1)

    ##########################
    #TODO write data to tiff #
    ##########################
    spueGeopoint = data2raster.GeoPoint(data['lon'],data['lat'],data['spue'])
    spueGeopoint.create_raster(filename="spue.tiff",output_format="GTiff")

    #########################
    # Get Cortad Data - SST #
    #########################
    ref_date=datetime.datetime(1980,12,31,12,0,0)

    dates = list()
    for data_date,data_time in zip(data['dates'],data['times']):
        data_datetime = data_date+data_time
        dates.append(datetime.datetime.strptime(data_datetime, '%y%m%d%H%M%S'))

    data_start = min(dates)
    data_end = max(dates)
    days = 60*60*24

    # calculate days from ref date to first sighting
    time_start = int(round((data_start - ref_date).total_seconds()/days))
    # calculate days from ref date to last sighting
    time_end = int(round((data_end - ref_date).total_seconds()/days))

    filledSST = getCortad.extractCORTADSST("North Sea",time_start,time_end)

    lonSST2D, latSST2D, lonSST, latSST = getCortad.extractCoRTADLongLat()

    print 'First sighting: ',data_start
    print 'Last sighting: ',data_end
    print 'Sighting Period: ', time_start, time_end, data_end-data_start

    #####################
    # Print Map Details #
    #####################

    lat0center,lon0center,mapWidth,mapHeight = centerMap(data['lat'],data['lon'],1.1)

    print 'Plotting Data From: '+data_file
    print "Map lat/lon center: ",lat0center,lon0center
    print "Map height/width: ",mapWidth,mapHeight


    ##############
    # Create Map #
    ##############
    # TODO write metadata to image file

    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111)

    # setup stereographic basemap.
    # lat_ts is latitude of true scale.
    # lon_0,lat_0 is central point.
    m = Basemap(width=mapWidth,height=mapHeight,
                resolution='l',projection='stere',\
                lat_0=lat0center,lon_0=lon0center)

    #m.shadedrelief()
    m.drawcoastlines(linewidth=0.2)
    m.fillcontinents(color='white', lake_color='aqua')

    # STT Plot
    drawSST(ax,m,lonSST2D,latSST2D,filledSST,0.25)

    ####################
    # plot effort grid #
    ####################
    #Setting the default values
    # TODO generalize this and re-add support for shapefiles, etc.
    xv, yv = m(data['lat'],data['lon'])
    data_values = data['spue']
    proj = None
    # TODO create ogr SpatialReference object for projection
    #proj= list('GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]')

    xSize= 200 # mapWidth
    ySize= 200 # mapHeight
    power=2.0
    smoothing=0.0
    driverName='GTiff'
    outFile = 'effort_interpolation'

    xMin,yMin = m(min(data['lat']),min(data['lon']))
    xMax,yMax = m(max(data['lat']),max(data['lon']))
    geotransform=[]
    geotransform.append(float(xMin))
    geotransform.append(float((xMax-xMin)/xSize))
    geotransform.append(float(0.0))
    geotransform.append(float(yMax))
    geotransform.append(float(0.0))
    geotransform.append(float((yMin-yMax)/ySize))
    print geotransform

    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,20.), labels=[0,0,0,1], fontsize=10)
    m.drawmapboundary(fill_color='aqua')

    plotSizedData(data['lat'],data['lon'],spue,'ro',4,20)
    #TODO get interpolation etc working
    #ti = 400.0
    #XI, YI = np.meshgrid(ti, ti)
    #ZI = invdistgis.invDist(xv,yv,data_values,geotransform,proj,xSize,ySize,power,smoothing,driverName,outFile)
    #plt.pcolor(XI, YI, ZI)

    # TODO cleanup leftovers from Trond script
    #plot data points
    #x, y = m(data['lon'],data['lat'])
    #m.scatter(x,y,2,marker='o',color='k')
    #plt.title("Example")
    #plotfile='SST_northsea_'+str(currentDate)+'.png'
    #print "Saving map to file %s"%(plotfile)
    #plt.savefig(plotfile)
    plt.show()
