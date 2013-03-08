#!/usr/bin/python
# -*- coding: utf-8 -*-
#Copyright: J Gomez-Dans <j.gomez-dans@geog.ucl.ac.uk>
#
#2009-08-18
#
#This code is licensed under the terms of the GNU General Public License (GPL),
#version 2 or above; See /usr/share/common-licenses/GPL , or see
#http://www.fsf.org/licensing/licenses/gpl.html
#
try:
    import osgeo.gdal as gdal
    import osgeo.ogr as ogr
    import osgeo.osr as osr
except:
    try:
        import gdal
        import ogr
        import osr
    except:
        import sys
        sys.stderr.write("GDAL needs to be installed. Exiting...\n")
        sys.exit(-1)

import numpy
import sys
import pdb
import os
import itertools
def SendMail(serverURL="MAILSERVER", sender='SENDER', to='', subject='',\
             text=''):
  import smtplib
  """
  Usage:
  mail('somemailserver.com', 'me@example.com', 
        'someone@example.com', 'test', 'This is a test')
  """
  headers = "From: %s\r\nTo: %s\r\nSubject: %s\r\n\r\n" %\
             (sender, to, subject)
  message = headers + text
  mailServer = smtplib.SMTP(serverURL)
  mailServer.sendmail(sender, to, message)
  mailServer.quit()

def ReadGrid ( file ):
        g = gdal.Open (file)
        grid = g.GetRasterBand(1).ReadAsArray().astype(numpy.int32)
        return grid

def CreateGDALDataSet ( fname ):
    from osgeo import gdal
    from osgeo import osr
    t_srs = osr.SpatialReference()
    t_srs.SetFromUserInput('WGS84')
    drv = gdal.GetDriverByName ("GTiff")
    ds = drv.Create ( fname, 360*6, 180*6, 16, gdal.GDT_Float32 )
    ds.SetProjection ( t_srs.ExportToWkt() )
    ds.SetGeoTransform ( [-180.0, 1./6., 0.0, 90.0, 0.0, -1./6.] )
    return ds

def WriteClassData ( mapa, fgrid, output_dir="/tmp/"):
  g = gdal.Open( fgrid )
  gridmap = g.GetRasterBand(1).ReadAsArray().astype(numpy.int32)
  for classification_scheme in mapa.iterkeys():
    canvas = numpy.zeros ((180*6,360*6,16),dtype=numpy.float32)
    fname = output_dir+'MODIS_classification_%s.tif'%classification_scheme
    ds = CreateGDALDataSet ( fname )
    for (grid_cell, histo) in mapa[classification_scheme].iteritems():
      total_pix = float(histo.sum())
      grid_where = gridmap==grid_cell
      for i in xrange(16):
    canvas[grid_where, i] = histo[i+1]/total_pix
      print
    for i in xrange(16):
      ds.GetRasterBand(i+1).WriteArray(canvas[:,:,i])
    ds = None


def RasteriseGrid ( grid_name, MODIS_fname):
    import os
    import tempfile
    (fd, grid_temp) = tempfile.mkstemp(suffix='.tif')
    os.close(fd)

    command_line = " gdal_translate -of GTiff -ot Int32 -scale 0 255 0 0.01"+\
                    " %s %s"%(MODIS_fname, grid_temp)
    os.system ( command_line )
    command_line = "gdal_rasterize -a grid_id -l word_1_6_grid_MODIS %s"+\
                    " %s"%(grid_name, grid_temp)
    os.system (  command_line )
    grid = ReadGrid ( grid_temp)
    #We create a dictionary with the grid ID as key. 
    #Each element has two arrays, the x pixel no and the y pixel no
    U = numpy.unique ( grid )
    R = dict ( itertools.izip ( U, [numpy.nonzero ( grid==U[i]) \
        for i in itertools.islice ( itertools.count(), U.shape[0])]))
    try:
        R.pop ( 0 ) # Do away with 0 valued pixels.
    except:
        pass
    os.remove ( grid_temp)
    return R



if __name__=="__main__":
    import glob
    import cPickle
    import sys, os, time
    schemas = ["IGBP","UMD","LAI","BGC","PFT"]
    files = glob.glob ("DIRECTORY_WITH_MODIS_MCD12Q1_FILES/MCD12Q1*.hdf")
    files.sort()
    spitfire_grid = "/MY_SHAPEFILE_WITH_GRId/word_1_6_grid_MODIS.shp"
    mapa = {}
    for lc_schema in [0]:#xrange(5):
        M={}
        for file in files:
            print schemas[lc_schema], file.split("/")[-1]
            #if (file.find("h17v05")>0) or ( file.find("h17v04")>0):
            fname =\
                'HDF4_EOS:EOS_GRID:"%s":MOD12Q1:Land_Cover_Type_1'%(file)
            print ">>> [%s] Starting rasterisation"%time.asctime()
            R = RasteriseGrid ( spitfire_grid, fname)
            print ">>> [%s] Finished rasterisation."%time.asctime()
            sys.stdout.flush()
            fname =\
            'HDF4_EOS:EOS_GRID:"%s":MOD12Q1:Land_Cover_Type_%1d'%(file, \
                    lc_schema+1)
            lc = ReadGrid ( fname )
            print ">>> [%s] Calculating histograms"%time.asctime()
            sys.stdout.flush()
            M_temp = dict(itertools.izip ( R.iterkeys(),\
                       [numpy.histogram ( lc[R[grid_cell]], \
                        bins=numpy.arange(0,256,1) )\
                        for grid_cell in R.iterkeys()]))
            print ">>> [%s] Finished histograms"%time.asctime()
            sys.stdout.flush()
            for (k,v) in M_temp.iteritems():
              M[k] = M.setdefault( k, v[0]) + v[0]
            print ">>> [%s] Finished dictionary update"%time.asctime()
            sys.stdout.flush()
        mapa[schemas[lc_schema]] = M
        WriteClassData ( mapa, "/GEOTIFF_WITH_RASTERISED_GRID",\
        output_dir="/OUTPUT_DIR/")
    mail_msg="""Hi!\n
    This is an automated email, to let you know that I have processed
    the MCD12Q1 data using the 1/6 degree grid. The results are to be
    found in </OUTPUT_DIR/>.

    This message was sent on %s
        """%time.asctime()

    SendMail(serverURL="SERVER", sender='SENDER', to='"NAME" EMAIL',\
         subject='Finished landcover aggregation!', text=mail_msg)

