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
        lons = nc_dataset.variables["lon"][:]
        lats = nc_dataset.variables["lat"][:]
        vals = nc_dataset.variables["CHL1_mean"][:][:]
#        date_list = np.zeros(len(nc_rows), dtype=(np.str_, 19))
        date_list[:] = str(datetime.datetime.strptime(file_date, '%Y%m%d'))

        lon_mesh, lat_mesh = np.meshgrid(lons,lats)

        nc_coord_ids = np.vstack([nc_rows,nc_cols])

        lons, lats = grid.isin_convert(coord=nc_coord_ids)
