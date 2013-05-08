geosetup
========

Python package for reprojecting, gridding, and writing geospatial data to formats supported by GDAL

**Code Credits:**
* Trond Kristiansen - plotting and CorTAD module
* more TBA

Modules
-------

* **writefile:** Given data and geo-points, writefile will write to file formats supported by the GDAL python, as well as a method for writing to ESRI's ASCII grid format (.asc).

* **gebco:** Extracts bathymetric data from gebco 1 degree and 30 arc-second global grids. Supports sub-setting given geographic bounds.

* **cortad:** Extracts SST data from weekly averaged cortad data. Modified from code written by Trond Kristiansen.

* **pathfinder:** Extracts SST from daily averaged pathfinder data. Supports subsetting by geographic and time bounds.

* **sightsurvey:** Imports data from textfile into structured numpy array for use with other modules.

* **makekml:** methods for creating Google kml files from dataset

* **interpolate:** Includes different interpolation methods for interpolating geo-spatial data. 

* **griddata:** Re-samples point/value data to new grid dimensions.

* **globcolour:** Extracts Chl-a data from daily averaged globcolour merged chlorophyll data. Supports subsetting by geographic and time bounds.


