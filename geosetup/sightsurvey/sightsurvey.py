import numpy as np
from netCDF4 import Dataset
import sys, os
from StringIO import StringIO
import datetime

#############
# Functions #
#############

def getData(data_file, skip_rows=0):

    with open(data_file) as fh:
        file_io = StringIO(fh.read().replace(',', '\t'))

    # Define names and data types for sighting data
    record_types = np.dtype([
                    ('vessel',str,1),           #00 - Vessel ID
                    ('dates',str,6),            #01 - Date
                    ('times',str,6),            #02 - Time (local?)
                    ('lat',float),              #03 - latitude dec
                    ('lon',float),              #04 - longitude dec -1
                    ('beafort',str,2),          #05 - beafort scale
                    ('weather',int),            #06 - weather code
                    ('visibility',int),         #07 - visibility code
                    ('effort_sec',float),       #08 - seconds on effort
                    ('effort_nmil',float),      #09 - n miles on effort
                    ('lat0',float),             #10 - lat of sight start
                    ('lon0',float),             #11 - lon of sight start
                    ('num_observers',int),      #12 - number of observers
                    ('species',str,6),          #13 - species codes
                    ('num_animals',int),        #14 - number observed
                    ('sighting',int),           #15 - Boolean sight code
                    ('rdist',float),            #16 - distance to sight
                    ('angle',float),            #17 - angel from ship
                    ('block',str,2),            #18 - cruise block
                    ('leg',int),                #19 - cruise leg code
                    ('observations',str,25),    #20 - cruise obs codes
                    ])

    # Import data to structured array
    data = np.genfromtxt(file_io, dtype = record_types, delimiter = '\t',
                         skip_header = skip_rows)

    # Correct longitude values
    data['lon'] = data['lon']*(-1)

    # Print a summary of geo data
    print '\nSighting Data Information'
    print '-------------------------------------------'
    print 'Data Path: ' + data_file
    print 'First sighting: ', min(data['dates'])
    print 'Last sighting: ', max(data['dates'])

    return data

    #TODO remove following
    # TODO calculate distance between start points and sighting point and
    # compare

    # Create array of unique survey block IDs
    #blocks = np.unique(data['block'])

    #g = pyproj.Geod(ellps='WGS84') # Use WGS84 ellipsoid
    #f_azimuth, b_azimuth, dist =
    #g.inv(data['lon'],data['lat'],data['lon0'],data['lat0'])

    #last_idx = 0
    #idx_pos = 0
    #minke_effort = np.zeros_like(minke_idx, dtype=float)
    #for idx in minke_idx:
    #    minke_effort[idx_pos] = dist[last_idx:(idx+1)].sum()
    #    last_idx = idx
    #    idx_pos = idx_pos+1

    # generate list of indexes where effort was greater than zero
    #effort_idx = np.where(data['effort_sec']*data['effort_nmil'] != 0)

    # spue = data['num_animals'][effort_idx]/(data['effort_sec'][effort_idx]*data['effort_nmil'][effort_idx])

    # append spue calculations to structured array dataset
    #data = np.lib.recfunctions.append_fields(data,'spue',data=spue)


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

    data_file = sys.argv[1]

    data, geobounds, timebounds = sightsurvey.getData(SIGHT_DATA, skip_rows=0)

    print 'data: ', data[:5]
    print 'geobounds: ', geobounds
    print 'timebounds: ', timebounds
