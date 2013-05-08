import pyproj
#import utm # faster than pyproj for utm, with auto UTM zone detection

'''http://all-geo.org/volcan01010/2012/11/change-coordinates-with-pyproj/'''

'''http://gis.stackexchange.com/questions/5683/basic-proj4-conversion-utm-grid-longlat-and-back/5684#5684'''

wgs84=pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
lcc=pyproj.Proj("+init=EPSG:3034") # Lambert Conformal Conical (LCC)
google_merc=pyproj.Proj("+init=EPSG:3857") # WGS84 Web Mercator (Auxillary Sphere; aka EPSG:900913)
osgb36=pyproj.Proj("+init=EPSG:27700") # UK Ordnance Survey, 1936 datum
utm26N=pyproj.Proj("+init=EPSG:32626") # UTM coords, zone 26N, WGS84 datum
utm27N=pyproj.Proj("+init=EPSG:32627") # UTM coords, zone 27N, WGS84 datum
utm28N=pyproj.Proj("+init=EPSG:32628") # ... you get the picture

# Define a projection with Proj4 notation, in this case an Icelandic grid
isn2004=pyproj.Proj("+proj=lcc +lat_1=64.25 +lat_2=65.75 +lat_0=65 +lon_0=-19 +x_0=1700000 +y_0=300000 +no_defs +a=6378137 +rf=298.257222101 +to_meter=1")

def transform_list(input_proj,output_proj,x,y):
    x2_list = list()
    y2_list = list()
    for point in zip(x,y):
        x2,y2 = pyproj.transform(input_proj,output_proj,point[0],point[1])
        x2_list.append(x2)
        y2_list.append(y2)
    return x2_list,y2_list

def transform_point(input_proj,output_proj,x,y):
    return pyproj.transform(input_proj,output_proj,x,y)

if __name__ == '__main__':
    # example coordinates
    lon = [-19.5, -19.7, -19.9,175.892]
    lat = [63.183, 63.583, 63.983, 45.232]
    p1 = wgs84
    p2 = google_merc

    print '\nInput x and y:\n',lon,lat

    print '\nOutput x and y:\n',transform_proj(p1,p2,lon,lat)

    x2,y2 = transform_list(p1,p2,lon,lat)
    print '\nOutput x and y to back to Input:\n',transform_proj(p2,p1,x2,y2),'\n'

#       print utm.from_latlon(lon[0],lat[0])
#       print utm.to_latlon(utm_x,utm_y,utm_zone_num,utm_zone_row)

# Expected output with isn2004:
# ([1674812.314071126, 1665231.4455360402, 1655933.043340466],
#  [97526.5926421243, 142220.30636996412, 186937.31348842022])
