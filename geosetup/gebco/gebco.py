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
    ''' Extract gebco bathymetry data.'''

    # Get Lat/Lon from first data file in list
    dataset = Dataset(file_path,'r')

    cols, rows = dataset.variables["dimension"]
    grid_w_deg, grid_h_deg = dataset.variables["spacing"]
    min_lon, max_lon = dataset.variables["x_range"]
    min_lat, max_lat = dataset.variables["y_range"]
    min_z, max_z = dataset.variables["z_range"]
    z = dataset.variables["z"][:5]

    dataset.close()

    print 'Gebco Bathymetric Data Summary:'
    print '------------------------------------'
    print 'cols: %i rows: %i' % (cols, rows)
    print 'grid width: %f grid height: %f' % (grid_w_deg, grid_h_deg)
    print 'min_lon: %f max_lon: %f' % (min_lon, max_lon)
    print 'min_lat: %f max_lat: %f' % (min_lat, max_lat)
    print 'min_z: %i max_z: %i' % (min_z, max_z)
    print 'First five z: ', z


def getGebcoData(file_path,min_lon,max_lon,min_lat,max_lat):
    '''Retrieve GEBCO depth over specified area

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

    # Calculate flattened index positions
    # Multiply number of elements in each row by row position, adding
    # the number of columns to get final element count
    # TODO Collecting all points between indices, even if center subset
#    idx_start = (min_lat_idx * cols) + min_lon_idx # (row1*z.shape[1])+col1
#    idx_end   = (max_lat_idx * cols) + max_lon_idx # (row2*z.shape[1])+col2
    idx = np.ravel_multi_index((np.mgrid[min_lat_idx:max_lat_idx,min_lon_idx:max_lon_idx].reshape(2,-1)), (rows,cols))
    # TODO remove
#    print 'idx_start: ', idx_start
#    print 'idx_end  : ', idx_end

    zi = 0
    z = np.zeros((max_lon_idx - min_lon_idx) * (max_lon_idx - min_lon_idx))
    for i in range(max_lat_idx - min_lon_idx):
        tmp = ((dataset.variables["z"][(i*cols):((i*cols)+cols)])[min_lon_idx:max_lon_idx])
        for j in tmp:
            z[zi] = j
            zi += 1

    # Extract data with flattened indexes
    #z = dataset.variables["z"][idx_start:idx_end]
    z = dataset.variables["z"][idx]
    dataset.close()

    # Create latitude and longitude arrays
    # TODO check if incorrect to offset by one
    lon_range = (max_lon_idx - min_lon_idx)+1
    lat_range = (max_lat_idx - min_lat_idx)+1
    lons = np.zeros(lon_range)
    lats = np.zeros(lat_range)

    for i in range(lon_range):
        lons[i] = ((i+min_lon_idx) - (360./cell_size/2))*cell_size

    for i in range(lat_range):
        lats[i] = ((i+min_lat_idx) - (180./cell_size/2))*cell_size

    # TODO remove
    print 'lon_range: ', lon_range
    print 'lat_range: ', lat_range
    print 'len_lons: ', len(lons)
    print 'len_lats: ', len(lats)
    print 'len_z:    ', len(z)

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
