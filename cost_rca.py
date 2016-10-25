# ---------------------------------------------------------------------------
# CREATE COST RCAS
# Created on: Jan 24 2006 by John Norman and David M. Theobald
# Last updated: March 18 2014 by Erin Peterson
# ---------------------------------------------------------------------------

# Import system modules
import arcgisscripting, sys, string, os, re, time #win32com.client
from time import *

# Create the Geoprocessor object
gp = arcgisscripting.create()

# Check out any necessary licenses
gp.CheckOutExtension("spatial")

# Script arguments...
#InDEM = sys.argv[1]
#InReach = sys.argv[2]
##InWaterBody = sys.argv[3]
#OutRCAShapeFile = sys.argv[4]
#WKSpace = sys.argv[5]

InDEM = "d:\\OneDrive\\work\\GIS\\Wen_streams\\elev_wen"
InReach = "d:\\OneDrive\\work\\research\\CHaMP\\GIS\\coverages\\Wenatchee\\WenSTHD_RV"
##InWaterBody = "c:\\temp\\wbgrid"
OutRCAShapeFile = "d:\\OneDrive\\work\\research\\CHaMP\\GIS\\coverages\\Wenatchee\\rca7"
WKSpace = "c:\\Temp"

try:
    # Get full path and Name for inputs
    InDEM = gp.Describe(InDEM).Path + "\\" + gp.Describe(InDEM).Name
    InReach = gp.Describe(InReach).Path + "\\" + gp.Describe(InReach).Name
    if InWaterBody <> "#":
            nWaterBody = gp.Describe(InWaterBody).Path + "\\" + gp.Describe(InWaterBody).Name

    gp.Workspace = WKSpace
    strWS = "RCA" + strftime("%Y%m%d%H%M%S", localtime())
    gp.CreateFolder ( gp.Workspace, strWS )
    gp.Workspace = gp.Workspace + "/" + strWS
    gp.AddMessage(" ")
    gp.AddMessage("Writing RCA temporary files to workspace: " + gp.Workspace )

    # set output extent of create HCA to the area displayed on the screen
    gp.overwriteoutput = 1


    # Perform zonalstats on water bodies and reaches to get the maximum reach OID value
    # within a waterbody.
    if InWaterBody <> "#":
        gp.AddMessage(" ")
        gp.AddMessage("Finding reach OIDs that match waterbodies....")
        gp.AddMessage("  ")
        if gp.Exists("burntemp"):
            gp.Delete("burntemp")
        gp.ZonalStatistics_sa(InWaterBody, "Value", InReach, "burntemp", "maximum", "DATA")
        InRaster = "burntemp"
    else:
        InRaster = InReach

    # Merge Burntemp and buffered reaches grids into seed10
    if gp.Exists("seed10"): gp.Delete("seed10")
    if InWaterBody <> "#":
        gp.AddMessage(" ")
        gp.AddMessage("Merging stream raster with new waterbodies to create seed raster")
        gp.AddMessage("  ")
        string = InReach + ";" + gp.Workspace + "\\" + "burntemp"
        
        bc = gp.GetRasterProperties_management(gp.Workspace + "\\burntemp", "BANDCOUNT")
        gp.MosaicToNewRaster(string, gp.Workspace, "seed10", "#", "#", "#", bc)
    else:
        gp.CopyRaster(InReach, gp.Workspace + "\\seed10")

    # Get the spatial extent of the c:/temp/seed10 raster for input into create contant raster
    desc = gp.Describe(InDEM)

    gp.Extent = "MAXOF"
    rasterExtent = desc.extent
    cellSize = desc.MeanCellHeight
    
    # Local variables...
    ann_raw = gp.Workspace + "\\ann_raw"
    bas10b = gp.Workspace + "\\bas10b"
    rca_ras = gp.Workspace + "\\rca_ras"
    Output_distance_raster = ""
    cost_back = gp.Workspace + "\\cost_back"
    flow_acc = gp.Workspace + "\\flow_acc"
    flow_lenb = gp.Workspace + "\\overlandist"
    slope_deg = gp.Workspace + "\\slope_deg"
    flow_weight = gp.Workspace + "\\flow_weight"
    ann_cost = gp.Workspace + "\\ann_cost"
    cost_surf = gp.Workspace + "\\cost_surf"
    flow_max = gp.Workspace + "\\flow_max"
    FlowDir_dem = gp.Workspace + "\\flowdir_dem_1"
    Instr_dir = gp.workspace + "\\strdir_gbl"
    Upstr_len = gp.Workspace + "\\flowaccum"
    RCAmaxFlow = gp.Workspace + "\\rcamax"
    RCAFlowdir = gp.Workspace + "\\rcadir"
    Instr_len  = gp.Workspace + "\\rcaflowlen"
    Output_drop_raster = ""

    # Process: annulus mean...
    gp.AddMessage("Finding ridgelines...")
    gp.SingleOutputMapAlgebra_sa("int((" + InDEM +" - focalmean(" + InDEM + ", annulus, 7, 12)) + 0.5) ", ann_raw, InDEM)

    # Process: Create Ann Weight...
    gp.AddMessage("Weighting ridgelines...")
    gp.SingleOutputMapAlgebra_sa("con(" + ann_raw + " > 0, pow(" + ann_raw + ", 1.90), 1)", ann_cost, ann_raw)

    if gp.Exists(ann_raw): gp.Delete(ann_raw)

    # Process: Slope...
    gp.AddMessage("Calculating slope (Degrees)...")
    gp.Slope_sa(InDEM, slope_deg, "DEGREE", "1")

    # Process: Flow Direction...
    gp.AddMessage("Calculating Flow Direction...")
    gp.FlowDirection_sa(InDEM, FlowDir_dem, "FORCE", Output_drop_raster)

    # Process: Flow Accumulation...
    gp.AddMessage("Calculating Flow Accumulation weights...")
    gp.FlowAccumulation_sa(FlowDir_dem, flow_acc, "")

    if gp.Exists(FlowDir_dem): gp.Delete(FlowDir_dem)

    # Process: Focal Statistics...
    gp.FocalStatistics_sa(flow_acc, flow_max, "Circle 6 CELL", "MAXIMUM", "DATA")

    if gp.Exists(flow_acc): gp.Delete(flow_acc)

    # Process: Single Output Map Algebra...
    gp.AddMessage("Calculating Hydro weights...")

    #gp.AddMessage(flow_weight + " = con(" + flow_max + " > 1, con(" + slope_deg + " > 1, 1.0 / (" + flow_max + " * tan(" + slope_deg + " DIV 57.2957)), 1), 1)")
    gp.MultiOutputMapAlgebra(flow_weight + " = con(" + flow_max + " > 1, con(" + slope_deg + " > 1, 1.0 / (" + flow_max + " * tan(" + slope_deg + " DIV 57.2957)), 1), 1)")

    # Process: Single Output Map Algebra (2)...
    gp.MultiOutputMapAlgebra(cost_surf + " =  con(" + ann_cost + " < 100000000, " + ann_cost + " * con(" + flow_weight + " >= 1, " + flow_weight + ", 0.1), 100000000)")

    # Process: Cost Allocation...
    gp.AddMessage("Generating Cost RCAs...")
    gp.CostAllocation_sa("seed10", cost_surf, rca_ras, "", "", "Value", Output_distance_raster, cost_back)

    if gp.Exists(slope_deg): gp.Delete(slope_deg)

    # Process: Overland flow and instream distance
    gp.AddMessage("Generating Overland Flow Raster...")
    str1 = "con([cost_back] <= 2, [cost_back], con([cost_back] == 3, 4, con([cost_back] == 4, 8, con([cost_back] == 5, 16, con([cost_back] == 6, 32, con([cost_back] == 7, 64, 128))))))"
    #gp.AddMessage("Generating Overland flow raster..." + str1)

    if gp.Exists("seed10"): gp.Delete("seed10")

    gp.SingleOutputMapAlgebra_sa(str1, bas10b, cost_back)

    # Process: Flow Length...
    gp.AddMessage("Calculating Flow Length Weights on Backlinks...")
    gp.FlowLength_sa(bas10b, flow_lenb, "DOWNSTREAM")

    # Process: Instream Flow length within a RCA...
    gp.FlowDirection_sa(flow_lenb, Instr_dir, "FORCE", Output_drop_raster)

    gp.FlowLength_sa(Instr_dir, Upstr_len, "UPSTREAM")

    gp.AddMessage("Building rca_ras VAT...")
    expression = "BuildVat " + rca_ras
    gp.MultiOutputMapAlgebra(expression)

    gp.ZonalStatistics(rca_ras, "Value", Upstr_len, RCAmaxFlow, "MAXIMUM", "DATA")

    str1 = RCAFlowdir + " = con(" + RCAmaxFlow + " == " + Upstr_len + ", 0, " + Instr_dir + ")"

    gp.MultiOutputMapAlgebra_sa(str1)
    gp.FlowLength_sa(RCAFlowdir, Instr_len, "DOWNSTREAM")


    # Process: Raster to Polygon...
    gp.AddMessage("Creating RCA Shapefile " + OutRCAShapeFile + " ...")
    gp.RasterToPolygon_conversion(rca_ras, "temp_rca.shp", "SIMPLIFY", "Value")
    gp.AddMessage("Dissolving Multi-part Polygons...")
    gp.Dissolve_management("temp_rca.shp", OutRCAShapeFile + ".shp", "GRIDCODE", "#", "MULTI_PART")

    # Process: Add Field...
    gp.AddField_management(OutRCAShapeFile + ".shp", "rca_id", "LONG", "16", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")

    # Process: Calculate Field...
    gp.AddMessage("Calculating RCA ID field rca_id for RCA shapefile " + OutRCAShapeFile + " ...")
    gp.CalculateField_management(OutRCAShapeFile + ".shp", "rca_id", "[GRIDCODE]")

    gp.AddMessage(" ")
    gp.AddMessage(" ")
    gp.AddWarning("Finished Create Cost RCAs Script")
    gp.AddMessage(" ")
    gp.AddMessage(" ")
except:
    gp.GetMessages()
