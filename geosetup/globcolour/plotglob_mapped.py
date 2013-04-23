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
Extract globcolour chlorophyll data from netCDFs
'''

#############
# Functions #
#############

def subset_geodata(max_lon, min_lon, max_lat, min_lat,lons, lats, values):
    '''subset_geodata  returns a subset of lat/lon/value data from defined
       defined bounds in decimal degrees'''
    subset_lons = lons[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    subset_lats = lats[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    subset_vals = values[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    return subset_lons, subset_lats, subset_vals

def find_nearest(array,value):
    '''
    Returns index position of element nearest to given value
    '''
    idx = (np.abs(array-value)).argmin()
    #return array[idx] # return value
    return idx # return index

def getMappedGlobcolour(data_dir, min_lon, min_lat, max_lon, max_lat, start_date, end_date):
    '''
    Extracts subset of globcolour data based on lat/lon bounds and dates in the
    iso format '2012-01-31'
    '''
    file_count = 0
    file_list = np.asarray(os.listdir(data_dir))
    file_dates = [re.split('[._]', nc_filename)[1] for filename in file_list]

    # Get Lat/Lon from first data file in list
    nc_dataset = Dataset(file_list[0],'r')
    lons = nc_dataset.variables["lon"][:]
    lats = nc_dataset.variables["lat"][:]

    # Get average of chlorophyl values from date range
    for data_file in file_list[date_idx]:
        current_file = os.path.join(data_dir,data_file)
        # import data as numpy array, then get columns
        nc_dataset = Dataset(data_file,'r')
        vals = nc_dataset.variables["CHL1_mean"][:][:]
        # create running sum of values
        if last_vals:
            vals = vals+last_vals
        last_vals = vals
        file_count += 1

    vals = vals/file_count #TODO simple average

    # date_list = np.zeros(len(nc_rows), dtype=(np.str_, 19))
    date_list[:] = str(datetime.datetime.strptime(file_date, '%Y%m%d'))

    # Mesh Lat/Lon the unravel to return lists
    #TODO make mesh optional (i.e. return grid or lists)
    lon_mesh, lat_mesh = np.meshgrid(lons,lats)
    nc_coord_ids = np.vstack([nc_rows,nc_cols])
    lons, lats = grid.isin_convert(coord=nc_coord_ids)

    # Print Globcolour Information
    print '\nChl-a Information'
    print '-------------------------------------------'
    print 'Start Date: ', start_date
    print 'End Date: ', end_date
    print 'Max Chl-a: ', np.amax(vals)
    print 'Min Chl-a: ', np.amin(vals)

    return lons, lats, chl_vals

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
    getMappedGlobcolour(sys.argv[1],0.,50.,30.,60.,'2007-05-01','2007-08-30')
