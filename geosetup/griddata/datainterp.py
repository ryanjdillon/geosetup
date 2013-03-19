import numpy as np
from scipy.interpolate import RectSphereBivariateSpline

def geointerp(lats,lons,data,grid_size_deg, mesh=False):
    '''We want to interpolate it to a global x-degree grid'''
    deg2rad = np.pi/180.
    new_lats = np.linspace(grid_size_deg, 180, 180/grid_size_deg)
    new_lons = np.linspace(grid_size_deg, 360, 360/grid_size_deg)
    new_lats_mesh, new_lons_mesh = np.meshgrid(new_lats*deg2rad, new_lons*deg2rad)

    '''We need to set up the interpolator object'''
   # lats = [float(lat) for lat in lats]
    print lats*2#*deg2rad
    lut = RectSphereBivariateSpline(lats*deg2rad, lons*deg2rad, data)

    '''Finally we interpolate the data. The RectSphereBivariateSpline
    object only takes 1-D arrays as input, therefore we need to do some reshaping.'''
    new_lats = new_lats_mesh.ravel()
    new_lons = new_lons_mesh.ravel()
    data_interp = lut.ev(new_lats,new_lons)

    if mesh == True:
        data_interp = data_interp.reshape((360/grid_size_deg,
                                           180/grid_size_deg)).T

    return new_lats/deg2rad, new_lons/deg2rad, data_interp


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    '''Suppose we have global data on a coarse grid'''
    lats = np.linspace(10, 170, 9) # in degrees
    lons = np.linspace(0, 350, 18) # in degrees
    data = np.dot(np.atleast_2d(90. - np.linspace(-80., 80., 18)).T,
                    np.atleast_2d(180. - np.abs(np.linspace(0., 350., 9)))).T

    # interpolate data to 1 degree grid
    new_lats, new_lons, data_interp = geointerp(lats,lons,data,10,mesh=False)
    print len(new_lats), len(data_interp)

    '''Looking at the original and the interpolated data,
    one can see that the interpolant reproduces the original data very well'''
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.imshow(data, interpolation='nearest')
    ax2 = fig.add_subplot(212)
    ax2.imshow(data_interp, interpolation='nearest')
    plt.show()

#       fig2 = plt.figure()
#       s = [3e9, 2e9, 1e9, 1e8]
#       for ii in xrange(len(s)):
#               lut = RectSphereBivariateSpline(lats, lons, data, s=s[ii])
#               data_interp = lut.ev(new_lats.ravel(),
#                               new_lons.ravel()).reshape((360, 180)).T
#               ax = fig2.add_subplot(2, 2, ii+1)
#               ax.imshow(data_interp, interpolation='nearest')
#               ax.set_title("s = %g" % s[ii])
#       plt.show()

## Plot gridded effort
## TODO sort the data by u
## http://docs.scipy.org/doc/numpy/reference/generated/numpy.sort.html
#data_lats = np.sort(data, order='lat')['lat']
#print data_lats
#data_lons = np.sort(data, order='lat')['lon']
#data_spue = np.sort(data, order='lat')['spue']
#lats_interp, lons_interp, data_interp = geointerp(data_lats,data_lons,data_spue,1)
#x2, y2 = map_object(data_lons,data_lats)
#levels=np.arange(2,18,0.5)
#
#map_object.contourf(x2,y2,data_interp,levels,cmap=plt.cm.jet,
#                                        extend='upper',alpha=myalpha)i
