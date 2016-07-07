year <- "2000"
setwd(paste0("D:/OneDrive/work/GIS/NARR/", year))
library(rgdal)
library(raster)
library(rNOMADS)
library(gdalUtils)
  
  grib <- readGDAL("narrmon-a_221_20000101_0000_000.grb")
  image(grib, attr=1)
  r2 <- raster("narrmon-a_221_20000101_0000_000.grb", layer="Temperature")
  r3 <- raster("narrmon-a_221_20000101_0000_000.grb", band=74)
  r4 <- raster("narrmon-a_221_20000101_0000_000.grb", band=269)
  r5 <- raster("narrmon-a_221_20000101_0000_000.grb", band=250)
  r6 <- stack("narrmon-a_221_20000101_0000_000.grb")
  
  gdalinfo(grib)
  
  gdal_setInstallation()
  valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
  if(require(raster) && require(rgdal) && valid_install)
  {
    gdal_translate("narrmon-a_221_20000101_0000_000.grb", "test_00", of = "NetCDF")
  }