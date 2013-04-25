import numpy as np
from netCDF4 import Dataset
import sys, os

#############
# Functions #
#############

def getGebco(data_dir, min_lon, min_lat, max_lon, max_lat):
    ''' Extract gebco bathymetry data.

    Uses the
    '''
    gebco_file = 'gebco_08.nc'
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
    z = dataset.variables["z"][:50]

    #print
    print 'cols: ',cols,' rows: ',rows
    print 'grid width: ', grid_w_deg
    print 'grid height: ', grid_h_deg
    print 'min_lon: ', min_lon
    print 'max_lon: ', max_lon
    print 'min_lat: ', min_lat
    print 'max_lat: ', max_lat
    print z_range
    print z


    dataset.close()

#    return lons, lats, depth

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

    getGebco(data_dir, 0., 40., 50., 65.)
