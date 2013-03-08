import numpy as np

'''
from: http://menugget.blogspot.no/2012/04/working-with-globcolour-data.html#more

This function is used converts ISIN grid information used by Globcolour to latitude
and longitude for a pefect sphere, as well as to construct associated polygons for use in mapping.

The raw Globcolour .nc files come with column and row pointers
as to the the grid's location. For 4.63 km resolution data, this
translates to 4320 latitudinal rows with varying number of
associated longitudinal columns depending on the latitudinal
circumference.

Input must be either a vector of grid numbers ["grd"] or a dataframe
with column and row identifiers ["coord", e.g. columns in coord$col and
rows in coord$row]

When the argument "polygon=FALSE" (Default), the function will output
a dataframe object containing grid information (gridnumber["output$grd"],
column["output$col"], row["output$row"], longitude["output$lon"], and
latitude[output$lat])

If the argument "polygon=TRUE", then the putput will be a list with
polygon shapes in a dataframe(longitudinal coordinates of corners
["[[i]]$x"], latitudinal coordinates of corners ["[[i]]$y"])
'''

def isin_convert(grid = None, coord = None, polygons = False):

    earth_radius = 6378.137 # TODO this was different in globcolour metadata
    Nlat = 4320 # Number of latitudinal bands, globcolour default 4km resolution
    pi = np.pi
    circum = 2*pi*earth_radius # circumference of Earth at equator

    lat_rows = np.arange(1,Nlat+1,1) # nparray of sequence 1 to Nlat

    # circumfrance at equator / lat bin width:
    # delta-radius, varies by how many parallels there are (Nlat)
    dr = (pi*earth_radius)/Nlat
    # angle increment between each parrallel:
    # delta-phi, from equator to pole
    dphi_lat = pi/Nlat

    phi = -(pi/2)+(lat_rows*dphi_lat)-(dphi_lat/2) # nparray of angles
    # calculate latitudal circumference of parallels:
    p = circum*np.cos(phi) # np array of circ in __ (2*Pi*r)*cos(phi)
    # number of lon rows, dimension = lat dimension at equ.
    Nlon = np.rint(p/dr) #round to nearest whole integer
    dlon = p/Nlon # circumfrance / number of meridians
    dphi_lon = (2*pi)/Nlon
    Ntot = sum(Nlon)
    lat = phi*(180/pi) # nparray of angles to radians

    if (grid is not None):
        # calculate coordinates
        cum_Nlon = np.cumsum(Nlon)
        # TODO the following doesn't work for array [1,2,3,4], matters?
        grid_rows = np.asarray([np.amax(np.where(cum_Nlon < point)) for point in grid])
        grid_cols = grid - cum_Nlon[grid_rows]
        # calculate longitude and latitude
        grid_lats = lat[grid_rows]
        Nlon_rows = Nlon[grid_rows]
        grid_lons = (360*(grid_cols-0.5)/Nlon_rows)-180

    if (coord is not None):
        #calculate coordinates
        grid_rows = coord[:][0]
        grid_cols = coord[:][1]
        print grid_rows
        print grid_cols
        # calcuate longitude and latitude
        grid_lats = lat[grid_rows]
        Nlon_rows = Nlon[grid_rows]
        grid_lons = (360*(grid_cols-0.5)/Nlon_rows)-180
        # calculate grid
        cum_Nlon = np.cumsum(Nlon)
        grid = cum_Nlon[grid_rows] + grid_cols # TODO check what this does
        cum_Nlon = np.cumsum(Nlon)
        return  grid_lats, grid_lons
    if (polygons is not None):
        Nlon_rows = Nlon[grid_rows]
        grid_width = 360/Nlon_rows
        grid_height = 180/Nlat
        # create list 1 to (number of elements in grid)
        polys = range(0,len(grid))
        xs = np.hstack([grid_lons[polys]-grid_width/2,
                        grid_lons[polys]-grid_width/2,
                        grid_lons[polys]+grid_width/2,
                        grid_lons[polys]+grid_width/2])

        ys = np.hstack([grid_lats[polys]-grid_height/2,
                        grid_lats[polys]+grid_height/2,
                        grid_lats[polys]+grid_height/2,
                        grid_lats[polys]-grid_height/2])
#TODO check that what's retu
#               return xs,ys #np.vstack([xs,ys])

    if (polygons is None):
        array = np.vstack([grid,grid_cols,grid_rows,grid_lons,grid_lats])

if __name__ == "__main__":
    grid = np.array([20,35,60,900])
    print isin_convert(grid)
