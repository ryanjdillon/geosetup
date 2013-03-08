# -*- coding: utf-8 -*-
#Copyright: J Gomez-Dans <j.gomez-dans@geog.ucl.ac.uk>
#
#2009-08-18
#
#This code is licensed under the terms of the GNU General Public License (GPL),
#version 2 or above; See /usr/share/common-licenses/GPL , or see
#http://www.fsf.org/licensing/licenses/gpl.html
#
from osgeo import gdal, osr
import numpy

g = gdal.Open("/tmp/world_grid.tif")
data = g.ReadAsArray()

a = numpy.ones ( (g.RasterXSize*g.RasterYSize), dtype=numpy.int32)
a = numpy.where ( data.flatten()==0, 0, a)
a = numpy.cumsum(a)
a = numpy.reshape (a,( g.RasterYSize, g.RasterXSize))
a = numpy.where ( data==0, 0, a)
a[900:,:] = 0

geo = [-180, 1./6, 0, 90., 0., -1./6]
nx = 360*6
ny = 180*6
driver = gdal.GetDriverByName ("GTiff")
dst_ds = driver.Create ( "/tmp/labelled_grid.tif",nx, ny, 1, gdal.GDT_Int32)
dst_ds.SetGeoTransform ( geo )
proj = osr.SpatialReference()
proj.ImportFromEPSG ( 4326)
dst_ds.SetProjection ( proj.ExportToWkt() )
dst_ds.GetRasterBand(1).WriteArray ( a)
dst_ds = None
