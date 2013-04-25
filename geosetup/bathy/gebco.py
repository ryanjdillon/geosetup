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
    idx = (coord/cell_size)+(deg_range/cell_size)

    # Round index to integer value up or down depending on boundary
    if limit_keyword == 'min':
        idx = math.floor(idx)
    elif limit_keyword == 'max':
        idx = math.ceil(idx)
    else :
        print '\nCoordinate keyword invalid. Exiting.\n'
        sys.exit()

    # Calculate decimal coordinate at rounded index position
    new_coord = (idx - (deg_range/cell_size))*cell_size

    return new_coord, idx

def summary(file_path, print_summary=False):
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

    if print_summary == True:
        print 'Gebco Bathymetric Data Summary:'
        print '------------------------------------'
        print 'cols: ', cols,             ' rows: ', rows
        print 'grid width: ', grid_w_deg, ' grid height: ', grid_h_deg
        print 'min_lon: ', min_lon,       ' max_lon: ', max_lon
        print 'min_lat: ', min_lat,       ' max_lat: ', max_lat
        print 'min_z: ', min_z,           ' max_z: ', max_z,
        print 'First five z: ', z


def getGebcoData(file_path,min_lon,max_lon,min_lat,max_lat):
    # Get Lat/Lon from first data file in list
    dataset = Dataset(file_path,'r')

    cell_size = dataset.variables["spacing"][0]

    min_lon, min_lon_idx = coord2idx(min_lon, 360, cell_size, 'min')
    max_lon, max_lon_idx = coord2idx(max_lon, 360, cell_size, 'max')
    min_lat, min_lat_idx = coord2idx(min_lat, 180, cell_size, 'min')
    max_lat, max_lat_idx = coord2idx(max_lat, 180, cell_size, 'max')

    print min_lon_idx
    print max_lon_idx
    print min_lat_idx
    print max_lat_idx

    z = dataset.variables["z"][min_lat_idx:max_lat_idx, min_lon_idx,max_lon_idx]
    dataset.close()

    # TODO Make the following better
    lons = np.array(len(z[0,:]))
    lats = np.array(len(z[:,0]))

    for i in range(len(lons)):
        lons[i] = min_lon + i(cell_size*(i+min_lon_idx))

    for i in range(len(lats)):
        lats[i] = min_lat + (cell_size*(i+min_lat_idx))

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

    lons, lats, z = getGebcoData(file_path,60,100,50,80)
