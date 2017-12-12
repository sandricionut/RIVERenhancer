
import os
import sys
import arcpy
from arcpy.sa import *
import multiprocessing
import numpy
import tempfile
import shutil

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
# Setup environment
arcpy.env.overwriteOutput = True

def xy2Point(x, y):
    return arcpy.Point(x, y)


def mCalculate(parameters):
    y = parameters[0][1]
    x = parameters[0][0]
    conditionalRaster = parameters[1]
    outFlowDirection = parameters[2]
    outputTextFile = parameters[3]
    counter = parameters[4]
    counter_lock = parameters[5]
    tempFolder = parameters[6]
    indexValue = None

    conditionalRasterArray = arcpy.RasterToNumPyArray(conditionalRaster)
    conditionalRasterUnique = numpy.unique(conditionalRasterArray)

    tempFolderProcess = os.path.join(tempFolder, "F_{}_{}".format(int(x), int(y)))

    if(os.path.exists(tempFolderProcess)):
        pass
    else:
        os.mkdir(tempFolderProcess)

    arcpy.env.scratchWorkspace = tempFolderProcess

    point = arcpy.Point(x, y)

    watershed = Watershed(outFlowDirection, arcpy.PointGeometry(point, arcpy.SpatialReference(31700)), None)
    watershedArray = arcpy.RasterToNumPyArray(watershed)
    watershedArrayUnique = numpy.unique(watershedArray)

    value0 = 0
    value1 = 0
    value3 = 0
    value5 = 0

    for valueWS in watershedArrayUnique:
        if(valueWS > 0):
            sumPixels = conditionalRasterArray[valueWS == watershedArray]
            for uniqueCR in conditionalRasterUnique:
                if(uniqueCR == 0):
                    value0 = int(numpy.sum(uniqueCR == sumPixels))
                if(uniqueCR == 1):
                    value1 = int(numpy.sum(uniqueCR == sumPixels))
                if(uniqueCR == 3):
                    value3 = int(numpy.sum(uniqueCR == sumPixels))
                if(uniqueCR == 5):
                    value5 = int(numpy.sum(uniqueCR == sumPixels))
        total = (value0 + value1 + value3 + value5)
        if(total > 0):
            indexValue = (value0 * 0.0 + value1 * 1.0 + value3 * 3.0 + value5 * 5.0) / total
            fText = open(outputTextFile, "a")
            fText.write("{}\t{}\t{}\n".format(point.X, point.Y, indexValue))
            fText.flush()
            fText.close()
        else:
            indexValue = -9999
    watershed = None
    shutil.rmtree(tempFolderProcess, ignore_errors =True)

    # counter_lock.acquire()
    # counter.value = counter.value + 1
    # counter_lock.release()
    # return indexValue


if __name__ == "__main__":

    try:
        tempFolder = os.path.join(tempfile.gettempdir(), "riverRenhancer")
        if(os.path.exists(tempFolder)):
            pass
        else:
            os.mkdir(tempFolder)

        arcpy.AddMessage("Temp folder is: " + tempFolder)

        # Input data
        dem = sys.argv[1] #arcpy.GetParameterAsText(0)
        conditionalRaster = sys.argv[2] #arcpy.GetParameterAsText(1)
        sampleSize = sys.argv[3] #arcpy.GetParameterAsText(2)
        samplePoints = sys.argv[4] #arcpy.GetParameterAsText(3)
        finalResult = sys.argv[5] #arcpy.GetParameterAsText(4)

        # Intermediate data - saved in temp folder
        streamRivers = os.path.join(tempFolder, "sr.tif")
        outFlowDir = os.path.join(tempFolder, "outFlowDir.tif")
        streamsRiversPoints = os.path.join(tempFolder, "streamRivers.shp")
        outputTextFile = os.path.join(tempFolder, "textRE.txt")

        # Setup environment
        arcpy.env.overwriteOutput = True
        arcpy.env.extent = conditionalRaster
        arcpy.env.snapRaster = conditionalRaster
        arcpy.env.scratchWorkspace = tempFolder

        _dem = Raster(dem)
        _nrCols = _dem.width
        _nrRows = _dem.height
        _cellSizeW = _dem.meanCellWidth
        _cellSizeH = _dem.meanCellHeight
        _minx = _dem.extent.XMin
        _maxy = _dem.extent.YMax
        _spatialReference = _dem.spatialReference

        # Execute FlowDirection
        outFlowDirection = FlowDirection(_dem, "FORCE", None)
        outFlowDirection.save(outFlowDir)
        # Execute FlowDirection
        outFlowAccumulation = FlowAccumulation(outFlowDirection, None, "INTEGER")
        # Execute Con using a map algebra expression instead of a where clause
        streamsRivers = Con(outFlowAccumulation >= 100, 1, 0)

        # Save streams and rivers pixels as points
        streamsRiversNull = arcpy.sa.SetNull(streamsRivers, streamsRivers, "VALUE = 0")
        arcpy.RasterToPoint_conversion(streamsRiversNull, streamsRiversPoints)

        streamsRivers.save(streamRivers)
        outFlowAccumulation = None
        outFlowDirection = None
        streamsRivers = None
        streamsRiversNull = None

        parameters = []
        outputTextFileTemp = open(outputTextFile, "w")
        outputTextFileTemp.close()

        # Define multiprocessing value and lock to keep track with model progress
        counter = multiprocessing.Value("i", 0)
        counter_lock = multiprocessing.Lock()

        for point in arcpy.da.SearchCursor(streamsRiversPoints, field_names=["SHAPE@XY"]):
            parameters.append([point[0], conditionalRaster, outFlowDir, outputTextFile, 0, 0, tempFolder])

        processes = multiprocessing.Pool(multiprocessing.cpu_count())
        result = processes.map(mCalculate, parameters)
        processes.join()
        print(result)

        # Process result
        arcpy.CreateFeatureclass_management(os.path.dirname(finalResult), os.path.basename(finalResult), geometry_type="POINT", spatial_reference=_spatialReference)
        arcpy.AddField_management(finalResult, field_name="indexValue", field_type="FLOAT", field_precision="8", field_scale="4")

        with arcpy.da.InsertCursor(finalResult, field_names=["SHAPE@XY", "indexValue"]) as rowInsert:
            for row in open(outputTextFile, "r").readlines():
                if(len(row.strip("\r\n").split("\t")) == 3):
                    print(row.strip("\r\n"))
                    x, y, index = row.strip().split("\t")
                    rowInsert.insertRow([arcpy.Point(x, y), index])


        # Check out the ArcGIS Spatial Analyst extension license
        arcpy.CheckInExtension("Spatial")

        arcpy.AddMessage("Finished multiprocessing")
    except arcpy.ExecuteError:
        # Geoprocessor threw an error
        arcpy.AddError(arcpy.GetMessages(2))
    except Exception as e:
        # Capture all other errors
        arcpy.AddError(str(e))