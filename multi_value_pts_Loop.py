 # ---------------------------------------------------------------------------
# PointToGridLoop.py
# Created on: Tue Jan 19 2010 10:59:23 AM
# Script to batch process point coverages of annual data to daily rasters
# ---------------------------------------------------------------------------
# Edit as needed for each directory/year/basin set

# Import arcpy module
import sys, string, os, arcgisscripting, arcpy

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Load required toolboxes...
gp.AddToolbox("C:\\Program Files (x86)\\ArcGIS\\Desktop10.2\\ArcToolbox\\Toolboxes\\Spatial Analyst Tools.tbx")

gp.workspace = "D:\\OneDrive\\work\\research\\CHaMP\\GIS\\LST\\LST_s1_2013\\"
outfolder = "D:\\OneDrive\\work\\research\\CHaMP\\GIS\\LST\\LST_s1_2013\\"
inPointFeatures = "D:\\OneDrive\\work\\research\\CHaMP\\GIS\\LST\\LST_s1_2013\\YF_13_LST.shp"

readFieldList = open("D:\\OneDrive\\work\\GIS\\8_day_1km_LST\\LST_s1_2013\\gridlist_num.txt", 'r')
FieldList = readFieldList.readlines()
readFieldList.close

for item in FieldList:
    item = item.replace("\n","")
    print item
    inraster = "lst_s1" + item 
    arcpy.gp.ExtractMultiValuesToPoints_sa(inPointFeatures, [[inraster, item]], "NONE")  
    
    
