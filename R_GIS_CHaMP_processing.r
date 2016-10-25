############################################################################################################
# This set of R scripts processes summer max temp model shapefiles for reprojection
# Created: 16 Oct 2016
# Used for all watersheds and all years by changing the basin names and year paths 
          

##############################################################################################

library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)
library(readxl)
library(lubridate)
library(zoo)
library(car)
library(gvlma)
library(ggplot2)
library(rgdal)
library(RColorBrewer)
library(classInt)

  basin <- "Ent"
  midBasin <- "Entiat"
  longBasin <- "Entiat"
  yrPath <- "12"
  yearPath <- "2012"
  dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/"
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)
  
  projname <- paste0(basin, "1km_poly_RCA_merge")
  
  setwd(paste0(dataPath, longBasin))
  
  proj_layer <- readOGR(dsn=".", layer = projname)
  
  
  netname <- paste0(basin, "_STHD_net")
  network <- readOGR(dsn=".", layer = netname)
  network <- spTransform(network, proj4string(proj_layer))
  
### change years and repeat this part #######
  
  setwd(paste0(mainPath, longBasin, "/", yearPath))
  
  preds <- read.dbf(paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf"))
  netmerge <- merge(network, preds, by.x='RCAID', by.y = 'RCAID')

  
# #### plot it to make sure it looks right ######### 
  
  
  fix3 <- classIntervals(netmerge@data[,4], n = 11, style = "fixed",fixedBreaks=c(6,8,10,12,14,15,16,18,20,22,24))
  fix3.colors <- findColours(fix3,pal=seis)
  plot(netmerge, col=fix3.colors, bg="black", fg="white")
  
# #################################################  
  
  setwd(paste0(mainPath, "shapes/"))
  layername <- paste0(basin, "_", yearPath, "_8D_summer_mx")
  writeOGR(netmerge, dsn=".", layer=layername, driver="ESRI Shapefile")
  
  
  