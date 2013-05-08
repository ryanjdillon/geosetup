import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
import matplotlib.pyplot as plt
from pylab import *
import mpl_util

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2012, 12, 3)
__modified__ = datetime.datetime(2012, 12, 3)
__version__  = "1.0"
__status__   = "Development, 03.12.2012"

def getCORTADtime():

    #base="http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version4/"
    base='/media/data/storage02/asf-fellowship/data/CoRTAD/version4/'
    file1="cortadv4_row00_col05.nc"
    filename1=base+file1
    cdf1=Dataset(filename1)

    print "Time: Extrating timedata from openDAP: %s"%(filename1)

    time=np.squeeze(cdf1.variables["time"][:])
    cdf1.close()
    return time

def openCoRTAD():
    # TODO Generalize to accept dates and geobounds, data dir
    """ Info on the different tiles used to identoify a region is found here:
    http://www.nodc.noaa.gov/SatelliteData/Cortad/TileMap.jpg"""
    #base='/media/data/storage02/asf-fellowship/data/CoRTAD/version4/'
    base= '/home/ryan/code/python/projects/asf/geosetup/data/cortad/'
#    base="http://data.nodc.noaa.gov/thredds/dodsC/cortad/Version4/"

    start_row = 00
    start_col = 05
    end_row = 01

    file1="cortadv4_row00_col05.nc"
    file2="cortadv4_row00_col06.nc"
    file3="cortadv4_row00_col07.nc"
    file4="cortadv4_row00_col08.nc"
    file5="cortadv4_row00_col09.nc"
    file6="cortadv4_row01_col05.nc"
    file7="cortadv4_row01_col06.nc"
    file8="cortadv4_row01_col07.nc"
    file9="cortadv4_row01_col08.nc"
    file10="cortadv4_row01_col09.nc"

    filename1=base+file1
    filename2=base+file2
    filename3=base+file3
    filename4=base+file4
    filename5=base+file5
    filename6=base+file6
    filename7=base+file7
    filename8=base+file8
    filename9=base+file9
    filename10=base+file10

    cdf1=Dataset(filename1)
    cdf2=Dataset(filename2)
    cdf3=Dataset(filename3)
    cdf4=Dataset(filename4)
    cdf5=Dataset(filename5)
    cdf6=Dataset(filename6)
    cdf7=Dataset(filename7)
    cdf8=Dataset(filename8)
    cdf9=Dataset(filename9)
    cdf10=Dataset(filename10)

    return cdf1,cdf2,cdf3,cdf4,cdf5,cdf6,cdf7,cdf8,cdf9,cdf10

def extractCoRTADLongLat():
    """Routine that extracts the longitude and latitudes for the
    combination of tiles. This is only necessary to do once so it is separated
    from the extraction of SST."""
#    cdf1,cdf2,cdf3,cdf4,cdf5,cdf6=openCoRTAD()
    cdf1,cdf2,cdf3,cdf4,cdf5,cdf6,cdf7,cdf8,cdf9,cdf10=openCoRTAD()

    longitude1=np.squeeze(cdf1.variables["lon"][:])
    latitude1=np.squeeze(cdf1.variables["lat"][:])

    longitude2=np.squeeze(cdf2.variables["lon"][:])
    latitude2=np.squeeze(cdf2.variables["lat"][:])

    longitude3=np.squeeze(cdf3.variables["lon"][:])
    latitude3=np.squeeze(cdf3.variables["lat"][:])

    longitude4=np.squeeze(cdf4.variables["lon"][:])
    latitude4=np.squeeze(cdf4.variables["lat"][:])

    longitude5=np.squeeze(cdf5.variables["lon"][:])
    latitude5=np.squeeze(cdf5.variables["lat"][:])

    longitude6=np.squeeze(cdf6.variables["lon"][:])
    latitude6=np.squeeze(cdf6.variables["lat"][:])

    longitude7=np.squeeze(cdf7.variables["lon"][:])
    latitude7=np.squeeze(cdf7.variables["lat"][:])

    longitude8=np.squeeze(cdf8.variables["lon"][:])
    latitude8=np.squeeze(cdf8.variables["lat"][:])

    longitude9=np.squeeze(cdf9.variables["lon"][:])
    latitude9=np.squeeze(cdf9.variables["lat"][:])

    longitude10=np.squeeze(cdf9.variables["lon"][:])
    latitude10=np.squeeze(cdf9.variables["lat"][:])

    cdf1.close();cdf2.close();cdf3.close();cdf4.close();cdf5.close();cdf6.close();cdf7.close();cdf8.close();cdf9.close();cdf10.close();

    longitude=concatenate((longitude1,longitude2,longitude3,longitude4,longitude5)) # 1-D array, axis irrelevant
    latitude=concatenate((latitude6,latitude1)) # 1-D array, axis irrelevant

    """ We have to flip this array so that we have increasing latitude
     values required by np.interp function. This means we also have to
     flip the input SST array"""
#    latitude=np.flipud(latitude) # flipping doesn not appear to be necessary
    lons,lats=np.meshgrid(longitude,latitude)

    print "Extracted longitude-latitude for CoRTAD region"
    print "Long min: %s Long max: %s"%(longitude.min(),longitude.max())
    print "Lat min: %s Lat max: %s"%(latitude.min(),latitude.max())
    print "------------------------------\n"
    return lons, lats, longitude, latitude

def extractCORTADSST(timestamp_start, timestamp_end, masked=True):
    """Routine that extracts the SST values for
    the specific tiles and time-period (t)"""
    cdf1,cdf2,cdf3,cdf4,cdf5,cdf6,cdf7,cdf8,cdf9,cdf10=openCoRTAD()

    # calculate days from ref date to first and last sighting
    ref_date=datetime.datetime(1980,12,31,12,0,0)
    days = 60.*60.*24. # sec*min*hr
    t1 = int(round((timestamp_start - ref_date).total_seconds()/days))
    t2 = int(round((timestamp_end - ref_date).total_seconds()/days))

    cortad_time=np.squeeze(cdf1.variables["time"][:])

    # use binary search to find index of nearest time value to data times
    idx1 = (np.abs(cortad_time-t1)).argmin()
    idx2 = (np.abs(cortad_time-t2)).argmin()

    # TODO make sure data is being averaged correctly
    filledSST1=np.average((cdf1.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST2=np.average((cdf2.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST3=np.average((cdf3.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST4=np.average((cdf4.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST5=np.average((cdf5.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST6=np.average((cdf6.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST7=np.average((cdf7.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST8=np.average((cdf8.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST9=np.average((cdf9.variables["FilledSST"][idx1:idx2,:,:]), axis=0)
    filledSST10=np.average((cdf10.variables["FilledSST"][idx1:idx2,:,:]), axis=0)

    offset=cdf1.variables["FilledSST"].__getattribute__('add_offset')
    cdf1.close();cdf2.close();cdf3.close();cdf4.close();cdf5.close();cdf6.close();cdf7.close();cdf8.close();cdf9.close();cdf10.close()

    filledMaskedSST1=filledSST1 - offset
    filledMaskedSST2=filledSST2 - offset
    filledMaskedSST3=filledSST3 - offset
    filledMaskedSST4=filledSST4 - offset
    filledMaskedSST5=filledSST5 - offset
    filledMaskedSST6=filledSST6 - offset
    filledMaskedSST7=filledSST7 - offset
    filledMaskedSST8=filledSST8 - offset
    filledMaskedSST9=filledSST9 - offset
    filledMaskedSST10=filledSST10 - offset

#    filledMaskedSST1=filledSST3*0
#    filledMaskedSST2=filledSST3*0
#    filledMaskedSST3=filledSST3*0
#    filledMaskedSST4=filledSST3*0
#    filledMaskedSST5=filledSST3*0
#    filledMaskedSST6=filledSST3*0

    """Now we have all the data in 4 different arrays that we need to concentate.
    First we add the horisontal tiles, and finally we stack the two horisontal ones on top
    of each other."""
    filledMaskedSST_lower=concatenate((filledMaskedSST1,filledMaskedSST2,filledMaskedSST3,filledMaskedSST4,filledMaskedSST5),axis=1)

    filledMaskedSST_upper=concatenate((filledMaskedSST6,filledMaskedSST7,filledMaskedSST8,filledMaskedSST9,filledMaskedSST10),axis=1)

    filledMaskedSST_all=concatenate((filledMaskedSST_upper,filledMaskedSST_lower),axis=0)

    """Flip the SST array to be consistent with order of latitude array"""
 #   filledMaskedSST_all=np.flipud(filledMaskedSST_all) # flipping doesn not appear to be necessary

    """ Scale and offset is autmoatically detected and edited by netcdf, but
    we need to mask the values that are not filled."""
    filledMaskedSST_final=ma.masked_less(filledMaskedSST_all,-2.)

    print "Min and max of SST: %s - %s"%(filledMaskedSST_final.min(),filledMaskedSST_final.max())
    print "------------------------------\n"

    return filledMaskedSST_final


if __name__ == "__main__":

    main()
