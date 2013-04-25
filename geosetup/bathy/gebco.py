import numpy as np
from netCDF4 import Dataset
import sys, os

#############
# Functions #
#############
def coord2idx(coord, deg_range, cell_size, limit_keyword, limits):
    '''Round decimal degrees or minutes to a defined limit'''

    # Convert decimal degress to components
    coord_degs = int(coord)
    coord_mins = (coord - int(coord))*60
    coord_secs = (coord_mins - int(coord_mins))*60

    # Set options depending on deg or seconds
    bin_size = cell_size/2. # (15*sec2deg)
    # TODO coord_val = 

    # check that coordinate is valid
    if (coord < -deg_range) or (coord > deg_range):
        print '\nBathymetric grid coordinate range incorrect. Exiting.\n'
        sys.exit()

    # Index offset giving value above or below, depending on keyword
    if limit_keyword == 'min':
        offset = 0
    elif limit_keyword == 'max':
        offset = 1
    else :
        print '\nCoordinate keyword invalid. Exiting.\n'
        sys.exit()

    def round2limit(coord_val, limits, offset):
        # Adjust coordinate value based on limits and offset
        # TODO generalize
        temp = 0.0
        for i in range(len(limits)-1):
            if (coord_val >= limits[i]):
                temp = limits[i + offset]
        coord_val = temp

    # Add components for adjusted coordinate
    new_coord = coord_degs+((coord_mins+(coord_secs/60.))/60.)

    # Correct for edge values overlapping meridian
    # TODO Doesn't make much sense with latitude
    if new_coord > deg_range:
        new_coord = new_coord - 360
    elif new_coord < -deg_range:
        new_coord = new_coord + 360

    # TODO finish
    if new_coord < 0:
        # reverse order from half array lenth
    elif new_coord >= 0:
        #add half array length

    # Divide by 0.5 cell width to get index position of coord
    idx_pos = (new_coord/(bin_size)) - 1 # subtract to get 0 position
    print 'index position:', idx_pos

    return new_coord, idx_pos

def getGebco(data_dir, min_lon, min_lat, max_lon, max_lat):
    ''' Extract gebco bathymetry data.

    Uses the
    '''
    #gebco_file = 'gebco_08.nc'
    # Array of second divisions in gebco file
    limits = np.array([0, 15, 30, 45, 60])

    gebco_file = 'gridone.nc'
    limits = np.array([0, 30, 60])

    # Get Lat/Lon from first data file in list
    dataset = Dataset(os.path.join(data_dir,gebco_file),'r')
    cols = dataset.variables["dimension"][0]
    rows = dataset.variables["dimension"][1]
    grid_w_deg = dataset.variables["spacing"][0]
    grid_h_deg = dataset.variables["spacing"][1]
    min_lon = dataset.variables["x_range"][0]
    max_lon = dataset.variables["x_range"][1]
    min_lat = dataset.variables["y_range"][0]
    max_lat = dataset.variables["y_range"][1]
    z_range = dataset.variables["z_range"][:]

    # Print Summary
    print 'cols: ',cols,' rows: ',rows
    print 'grid width: ', grid_w_deg
    print 'grid height: ', grid_h_deg
    print 'min_lon: ', min_lon
    print 'max_lon: ', max_lon
    print 'min_lat: ', min_lat
    print 'max_lat: ', max_lat
    print z_range
    print 'lenght z: ', len(dataset.variables["z"][:])

    dataset.close()


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

    limits = np.array([0, 15, 30, 45, 60])

    test_coord = -179.45869
    print coord2idx(test_coord, 'max', limits)
    print coord2idx(75.75, 'min', limits)
    print coord2idx(180., 'max', limits)
    print coord2idx(180., 'min', limits)
    getGebco(data_dir, 0., 40., 50., 65.)
