import sys, os
from netCDF4 import Dataset
import numpy as np
import isingrid as grid
import re
import csv
import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
'''
extract globcolour chlorophyll data from netCDFs
__author__: Ryan J. Dillon
'''

#############
# Functions #
#############

def subset_coord_data(max_lon, min_lon, max_lat, min_lat,lons, lats, values):
    '''subset_coord_data  returns a subset of lat/lon/value data from defined
       defined bounds in decimal degrees'''
    subset_lons = lons[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    subset_lats = lats[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    subset_vals = values[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    return subset_lons, subset_lats, subset_vals

################
# Main Program #
################

if __name__ == '__main__':

    #####################
    # commandline usage #
    #####################
    if len(sys.argv) < 2:
        print >>sys.stderr,'Usage:',sys.argv[0],'<data directory>\n'
        sys.exit(1)

#       nc_file = '/home/ryan/Desktop/asf-fellowship/code/globcolour/L3b_20070101__GLOB_4_GSM-MERMODSWF_CHL1_DAY_00.nc'

#    print sys.argv[1:]

    for datafile in sys.argv[1:]: # must provide dir with wildcard: ./*nc
        print datafile
        nc_filepath, nc_filename = os.path.split(datafile)
        file_date = re.split('[._]', nc_filename)[1]

        # import data as numpy array, then get colums
        #nc_dataset = Dataset(nc_file,'r')
        nc_dataset = Dataset(datafile,'r')
        nc_rows = nc_dataset.variables["row"][:]
        nc_cols = nc_dataset.variables["col"][:]
        nc_vals = nc_dataset.variables["CHL1_mean"][:]
        date_list = np.zeros(len(nc_rows), dtype=(np.str_, 19))
        date_list[:] = str(datetime.datetime.strptime(file_date, '%Y%m%d'))

        nc_coord_ids = np.vstack([nc_rows,nc_cols])

        lons, lats = grid.isin_convert(coord=nc_coord_ids)
        #print lon.shape,lat.shape,nc_vals.shape

        #print 'lon\n',np.amax(lon), np.amin(lon), lon.shape
        #print 'lat\n',np.amax(lat), np.amin(lat), lat.shape
        #print 'ncvals\n',np.amax(nc_vals), np.amin(nc_vals), np.size(nc_vals)
        #print file_date
        #print date_list
        #print np.swapaxes(np.vstack([lon,lat,nc_vals]),0,1)

        #nc_all = np.vstack([file_date,nc_coord_ids,nc_vals])
        #grid_coords = isin.isin_convert(nparray_coords)
        #grid_chla =

        # Write CSV
        #coordFile = open('./test-coords.csv', 'wb')
        #coordWriter = csv.writer(coordFile)
        #coordWriter.writerows(np.swapaxes(np.vstack([lon,lat,nc_vals]),0,1))
        #print np.swapaxes(np.vstack([lon,lat]),0,1)
        #coordFile.close()

        #########################
        # Create and fill array #
        #########################
        #subset_coord_data(max_lon, min_lon, max_lat, min_lat,lons,lats,values)
        lons, lats, vals = subset_coord_data(30,-15,65,40,lons,lats,nc_vals)

        ##################
        # Reproject Data #
        ##################

        # LatLon with WGS84 datum used by GPS units and Google Earth
        wgs84 = pyproj.Proj("+init=EPSG:4326")
        # Lambert Conformal Conical (LCC)
        lcc = pyproj.Proj("+init=EPSG:3034")
        # WGS84 Web Mercator (Auxillary Sphere; aka EPSG:900913)
        web_mercator = pyproj.Proj("+init=EPSG:3857")

        ##############
        # Create Map #
        ##############
        m = Basemap(width=12000000,height=8000000,
                    resolution='l',projection='stere',
                    lat_0=60,lon_0=0)

        m.drawcoastlines(linewidth=0.2)
        m.fillcontinents(color='white', lake_color='aqua')

        x_lon, y_lat = m(lons,lats)
        x_lon, y_lat = np.asarray(x_lon), np.asarray(y_lat)

        N, M = len(np.unique(x_lon)), len(np.unique(y_lat))

        result = np.empty((N,M))
        result.fill(np.nan)
        result.flat[:len(vals)] = vals

        levels=np.arange(2,18,0.5)
        CS2 = m.contourf(x_lon,y_lat,result,levels,cmap=plt.cm.jet,extend='upper')
        plt.show()
