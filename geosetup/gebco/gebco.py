import numpy as np
from netCDF4 import Dataset
import sys, os
import math

#############
# Functions #
#############

def coord2idx(coord, deg_range, cell_size, limit_keyword):
    '''Round decimal degrees or minutes to a defined limit'''

    # check that coordinate is valid
    if (coord < -deg_range) or (coord > deg_range):
        print '\nBathymetric grid coordinate range incorrect. Exiting.\n'
        sys.exit()

    # Calculate index position from provided coordinate
    idx = (coord/cell_size)+(deg_range/cell_size/2)

    # Round index to integer value up or down depending on boundary
    if limit_keyword == 'min':
        idx = int(math.floor(idx))
    elif limit_keyword == 'max':
        idx = int(math.ceil(idx))
    else :
        print '\nCoordinate keyword invalid. Exiting.\n'
        sys.exit()

    # Calculate decimal coordinate at rounded index position
    new_coord = (idx - (deg_range/cell_size/2))*cell_size

    return new_coord, idx

def summary(file_path):
    ''' Extract gebco bathymetry data summary.'''

    # Get Lat/Lon from first data file in list
    dataset = Dataset(file_path,'r')

    cols, rows = dataset.variables["dimension"]
    grid_w_deg, grid_h_deg = dataset.variables["spacing"]
    min_lon, max_lon = dataset.variables["x_range"]
    min_lat, max_lat = dataset.variables["y_range"]
    min_z, max_z = dataset.variables["z_range"]
    z = dataset.variables["z"][:5]

    dataset.close()

    print '\nGebco Bathymetric Data Summary:'
    print '------------------------------------'
    print 'cols: %i rows: %i' % (cols, rows)
    print 'grid width: %06.5f grid height: %06.5f' % (grid_w_deg, grid_h_deg)
    print 'min_lon: %5.1f max_lon: %5.1f' % (min_lon, max_lon)
    print 'min_lat: %5.1f max_lat: %5.1f' % (min_lat, max_lat)
    print 'min_z: %i max_z: %i' % (min_z, max_z)
    print 'First five z: ', z


def getGebcoData(file_path,min_lon,max_lon,min_lat,max_lat):
    '''Extract gebco bathymetric data from geographic bounds

       cell_size: decimal degree x & y dimension of grid cells
       The depth data is a 1-D array. Given its size, it is faster
       to use fancy indexing to extract the geographical subsection.
    '''
    dataset = Dataset(file_path,'r')
    cell_size = dataset.variables["spacing"][0]
    cols, rows = dataset.variables["dimension"]

    # Set lon (column) index and lat (row) index, retrieve adjusted coords 
    min_lon, min_lon_idx = coord2idx(min_lon, 360, cell_size, 'min')
    max_lon, max_lon_idx = coord2idx(max_lon, 360, cell_size, 'max')
    min_lat, min_lat_idx = coord2idx(min_lat, 180, cell_size, 'min')
    max_lat, max_lat_idx = coord2idx(max_lat, 180, cell_size, 'max')

    # TODO remove
    print 'min_lon_idx', min_lon_idx, 'min_lon', min_lon
    print 'max_lon_idx', max_lon_idx, 'max_lon', max_lon
    print 'min_lat_idx', min_lat_idx, 'min_lat', min_lat
    print 'max_lat_idx', max_lat_idx, 'max_lat', max_lat

    lon_range = max_lon_idx - min_lon_idx
    lat_range = max_lat_idx - min_lat_idx
    data_range = lon_range * lat_range

    zi = 0
    # Create zero array with the appropriate length for the data subset
    z = np.zeros(data_range)
    print z.shape
    # Process number of rows for which data is being extracted
    for i in range((max_lat_idx - min_lat_idx)):
        # Pull row, then desired elements of that row into buffer
        tmp = (dataset.variables["z"][(i*cols):((i*cols)+cols)])[min_lon_idx:max_lon_idx]
        # Add each item in buffer sequentially to data array
        for j in tmp:
            z[zi] = j
            # Keep a count of what index position the next data point goes to
            zi += 1

    dataset.close()

    # Create latitude and longitude arrays
    # TODO check if incorrect to offset by one
#    lons = np.empty(data_range).fill(min_lon_idx)
#    lats = np.zeros(data_range).fill(min_lat_idx)
    lon_const = 360./cell_size/2
    lat_const = 180./cell_size/2
    lons = (np.asarray(range(lon_range)) + min_lon_idx - lon_const)*cell_size
    lats = (np.asarray(range(lat_range)) + min_lat_idx - lat_const)*cell_size

#    for i in range(data_range):
#        lons = ((i+min_lon_idx) - (360./cell_size/2))*cell_size

#    for i in range(data_range):
#        lats[i] = ((i+min_lat_idx) - (180./cell_size/2))*cell_size

    # TODO remove
    print 'len_lons: ', lons.shape
    print 'len_lats: ', lats.shape
    print 'len_z:    ', z.shape

    return lons, lats, z

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
    data_file = 'gridone.nc'

    file_path = os.path.join(data_dir,data_file)

    summary(file_path)

                  # getGebcoData(file_path,min_lon,max_lon,min_lat,max_lat):
    lons, lats, z = getGebcoData(file_path,-180,0,0,90)

    print 'lons: ', lons[:]
    print 'lats: ', lats[:]
    print 'z: ', z[:]
