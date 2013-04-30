import sys, os
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
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

def printstuff(vals,vals_sum,file_count):
    # create running sum of values
    print file_count
    print '----------------'
    print 'sum', ma.sum(vals)
    print 'total sum', ma.sum(vals_sum)
    print '----------------'
    print 'max', ma.amax(vals)
    print 'min', ma.amin(vals)
    print '\n----------------'

def getMappedGlobcolour(data_dir, min_lon, min_lat, max_lon, max_lat, start_date, end_date):
    '''
    Extracts subset of globcolour data based on lat/lon bounds and dates in the
    iso format '2012-01-31'
    '''
    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d')
    file_list = np.asarray(os.listdir(data_dir))
    file_dates = np.asarray([datetime.datetime.strptime(re.split('[._]', filename)[1], '%Y%m%d')
                             for filename in file_list])
    data_files = file_list[(file_dates >= start_date) & (file_dates <= end_date)]

    file_count = 0

    # Get Lat/Lon from first data file in list
    dataset = Dataset(os.path.join(data_dir,file_list[0]),'r')
    lons = dataset.variables["lon"][:]
    lats = dataset.variables["lat"][:]
    vals_sum = np.zeros_like(dataset.variables["CHL1_mean"][:][:])
    mask_sum = np.empty(np.shape(vals_sum),dtype=bool)
    dataset.close()

    # Create cumulative mask before averaging data
    for data_file in data_files:
        current_file = os.path.join(data_dir,data_file)
        dataset = Dataset(current_file,'r') # by default numpy masked array
        mask = (dataset.variables["CHL1_mean"][:][:]).mask
        mask_sum =  mask + mask_sum

    vals_sum.mask = mask_sum

    # Get average of chlorophyl values from date range
    for data_file in data_files:
        current_file = os.path.join(data_dir,data_file)
        dataset = Dataset(current_file,'r') # by default numpy masked array
        vals = dataset.variables["CHL1_mean"][:][:]
        vals.mask = mask_sum
        vals_sum += vals
        #TODO remove
        #printstuff(vals,vals_sum,file_count)
        file_count += 1
        dataset.close()

    vals_mean = vals/file_count #TODO simple average

    # Mesh Lat/Lon the unravel to return lists
    #TODO make mesh optional (i.e. return grid or lists)
    lons_mesh, lats_mesh = np.meshgrid(lons,lats)
    lons = np.ravel(lons_mesh)
    lats = np.ravel(lats_mesh)

    # Print Globcolour Information
    print '\nChl-a Information'
    print '-------------------------------------------'
    print 'Start Date: ', start_date
    print 'End Date: ', end_date
    print 'Max Chl-a: ', np.amax(vals_mean)
    print 'Min Chl-a: ', np.amin(vals_mean)

    return lons, lats, vals_mean

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

    data_dir = sys.argv[1]

    lons, lats, vals_mean = getMappedGlobcolour(data_dir,0.,50.,30.,60.,'2007-08-01','2007-08-30')
    print lons
