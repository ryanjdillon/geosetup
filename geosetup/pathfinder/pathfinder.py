import sys, os
import math
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import re
import csv
import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
'''
Extract pathfinder sst data from netCDFs

this is a copy of globcolour.py and perhaps a generalized script for handling
both may be more desirable.
'''

#############
# Functions #
#############

def subset_geodata(max_lon, min_lon, max_lat, min_lat,lons, lats, values):
    '''Returns a subset of lat/lon/value data 

       from defined defined bounds in decimal degrees'''

    subset_lons = lons[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    subset_lats = lats[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    subset_vals = values[(lons<max_lon)&(lons>min_lon)&(lats<max_lat)&(lats>min_lat)]
    return subset_lons, subset_lats, subset_vals

def find_nearest(array,value):
    '''Returns index position of element nearest to given value'''

    idx = (np.abs(array-value)).argmin()
    #return array[idx] # return value
    return idx # return index

def printstuff(vals,vals_sum,file_count):
    # create running sum of values
    print file_count
    print '----------------'
    print 'sum', ma.sum(vals)
    print 'total sum', ma.sum(vals_sum)
    print 'perc masked:', ma.count_masked(vals_sum)/vals_sum.size*100.,'%'
    print '----------------'
    print 'max', ma.amax(vals)
    print 'min', ma.amin(vals)
    print '\n----------------'

def getnetcdfdata(data_dir, nc_var_name, min_lon, max_lon, min_lat, max_lat, data_time_start, data_time_end):
    '''Extracts subset of data from netCDF file 

    based on lat/lon bounds and dates in the iso format '2012-01-31' '''

    # Create list of data files to process given date-range
    file_list = np.asarray(os.listdir(data_dir))
    file_dates = np.asarray([datetime.datetime.strptime(re.split('-', filename)[0], '%Y%m%d%H%M%S') for filename in file_list])
    data_files = np.sort(file_list[(file_dates >= data_time_start)&(file_dates <= data_time_end)])

    # Get Lat/Lon from first data file in list
    dataset = Dataset(os.path.join(data_dir,file_list[0]),'r')
    lons = dataset.variables["lon"][:]
    lats = dataset.variables["lat"][:]

    # Create indexes where lat/lons are between bounds 
    lons_idx = np.where((lons > math.floor(min_lon))&(lons < math.ceil(max_lon)))[0]
    lats_idx = np.where((lats > math.floor(min_lat))&(lats < math.ceil(max_lat)))[0]
    x_min = lons_idx.min()
    x_max = lons_idx.max()
    y_min = lats_idx.min()
    y_max = lats_idx.max()

    # Create arrays for performing averaging of files
    vals_sum = np.zeros((y_max-y_min + 1, x_max-x_min + 1))
    vals_sum = ma.masked_where(vals_sum < 0 , vals_sum)
    mask_sum = np.empty((y_max-y_min + 1, x_max-x_min + 1))
    dataset.close()

    #TODO correct masking/interpolation.
    # Create cumulative mask for averaging data in provided date range
#    for data_file in data_files:
#        current_file = os.path.join(data_dir,data_file)
#        dataset = Dataset(current_file,'r') # by default numpy masked array
#        mask = (dataset.variables[nc_var_name][0 , y_min:y_max, x_min:x_max]).mask
#        mask_sum =  mask + mask_sum
#    vals_sum.mask = mask_sum

    # Get average of chlorophyl values from date range
    file_count = 0
    for data_file in data_files:
        current_file = os.path.join(data_dir,data_file)
        dataset = Dataset(current_file,'r') # by default numpy masked array
        # TODO filter this by quality flags
        vals = np.copy(dataset.variables[nc_var_name][0 ,y_min:y_max + 1, x_min:x_max + 1])
        # TODO correct so masking/interpolating instead of zeros
        #vals = ma.masked_where(vals < 0, vals)
        #ma.set_fill_value(vals, -999)
        vals = vals.clip(0) #TODO remove once above is corrected
        vals_sum += vals
        file_count += 1
        dataset.close()
    vals_mean = vals_sum/file_count #TODO simple average

    # Mesh Lat/Lon the unravel to return lists
    #TODO make mesh optional (i.e. return grid or lists)
    lons_mesh, lats_mesh = np.meshgrid(lons[lons_idx],lats[lats_idx])
    lons = np.ravel(lons_mesh)
    lats = np.ravel(lats_mesh)
    vals_mean = np.ravel(vals_mean)

    # Print Globcolour Information
    print '\n'+nc_var_name+'SST Information'
    print '-------------------------------------------'
    print 'Start Date:', data_time_start
    print 'End Date:', data_time_end
    print 'Max '+nc_var_name+':', np.amax(vals_mean)
    print 'Min '+nc_var_name+':', np.amin(vals_mean)

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
    start_date = '2007-08-01'
    end_date = '2007-08-30'
    start_date = datetime.datetime.strptime(start_date, '%Y-%m-%d')
    end_date = datetime.datetime.strptime(end_date, '%Y-%m-%d')

    lons, lats, vals_mean = getnetcdfdata(data_dir,'sea_surface_temperature',0.,50.,30.,60.,
                                          start_date, end_date)
    print lons[0]-lons[1]
    print lons[1]-lons[2]
    print lons[100]-lons[101]
