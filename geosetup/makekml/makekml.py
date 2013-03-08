import sys 
import csv 
import simplekml

# commandline usage
print '\n'
if len(sys.argv) < 3:
        print >>sys.stderr,'Usage:',sys.argv[0],'<datafile> <#rows to skip>'
        sys.exit(1)

dataFile = sys.argv[1]
dataStream = open(dataFile, 'rb')
dataReader = csv.reader(dataStream,delimiter='\t')
skip_rows = int(sys.argv[2])

#skipRows(dataReader, numRows) #should make this option, and default to none

dataValues = []
dataLat = []
dataLon = []

print 'Creating KML From: '+dataFile

[dataReader.next() for row in range(skip_rows)]
for row in dataReader:
#        dataValues.append(float(row[21]))
        dataLat.append(float(row[4]))
        dataLon.append(float(row[6]))
	
print 'lat '+str(max(dataLat)),'lon '+str(max(dataLon))
print 'lat '+str(min(dataLat)),'lon '+str(min(dataLon))
# Create an instance of Kml
kml = simplekml.Kml(open=1)

# Create a point named "The World" attached to the KML document with its coordinate at 0 degrees latitude and longitude.
# All the point's properties are given when it is constructed.
#single_point = kml.newpoint(name="The World", coords=[(0.0,0.0)])

# Create a point for each lat/lon pair. Point properties assigned after point created
for lat, lon in zip(dataLon,dataLat):
	pnt = kml.newpoint()
	pnt.coords = [(lat, lon)]
# Save the KML
kml.save(sys.argv[1]+'.kml')
