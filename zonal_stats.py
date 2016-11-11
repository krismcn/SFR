# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# zonal_stats.py
# Created on: 2016-06-29 13:13:42.00000
#   (generated by ArcGIS/ModelBuilder)
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import sys, string, os, arcgisscripting, arcpy

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Load required toolboxes...
gp.AddToolbox("C:\\Program Files (x86)\\ArcGIS\\Desktop10.2\\ArcToolbox\\Toolboxes\\Spatial Analyst Tools.tbx")

gp.workspace = "D:\\OneDrive\work\\research\\GIS\\NARR\\2000\\Temp"
outfolder = "D:\\OneDrive\\work\\research\\GIS\\NARR\2000\\Temp"

readFieldList = open("D:\\OneDrive\\work\\GIS\\NARR\\2000\\Temp\\gridlist.txt", 'r')
FieldList = readFieldList.readlines()
readFieldList.close

# Local variables:
Outline = "D:\\Dropbox\\work\\research\\CHaMP\\GIS\\coverages\\Wenatchee\\Wen_basin.shp"
startYear = 0101
# Process: Zonal Statistics as Table
for item in FieldList:
    item = item.replace("\n","")
    print item
    outFile = "Wen_Temp_" + str(startYear)
    outFile2 = outFile + ".dbf"
    arcpy.gp.ZonalStatisticsAsTable_sa(Outline, "Id", item, outFile, "DATA", "ALL")
    startYear = startYear + 1
    

                    