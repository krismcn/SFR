############################################
# R code to read in GIS data and summarize over specified time periods
# Created 17 Nov 2016
############################################

library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)
library(readxl)
library(zoo)
library(car)
library(gvlma)
library(ggplot2)
library(rgdal)
library(RColorBrewer)
library(classInt)
library(spatialEco)
library(maptools)
library (rgeos)


dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"

setwd(paste0(mainPath, "Summary_metrics/"))

  mikesPts <- read.dbf("Mikespts.dbf")
  mikesPts$SiteName <- as.character(mikesPts$SiteName)
  
setwd(paste0(dataPath, "All_CHaMP/"))

  ptsName <- "IC_2GIS2_SitesInCHamP_pj"
  pts <- readOGR(dsn=".", layer=ptsName)

  ptIDs <- pts@data[,c(13,14,45, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81)]
  ptIDs$Site_ID <- as.character(ptIDs$Site_ID)
  ptIDs[ptIDs == 0] <- NA
  
  basins <- as.character(unique(pts@data$WatershedN))

# ################
# This section repeats for each basin/year
# ################
basin <- "YF"
longBasin <- "YankeeFork"
yrPath <- 13
yearPath <- 2013 

setwd(paste0(mainPath, yearPath, "_temp_CHaMP/"))
  
  #modelName <- paste0("predt", yearPath, "_", basin, "_8D_mn")
  modelName <- paste0(basin, "_", yearPath, "_8D_mn")
  tempEst <- read.dbf(paste0(modelName, ".dbf"))
  colnames(tempEst)[48] <- "RCAID"  
  tempEst[1:5,29:33]
  tempEst[tempEst == -9999] <- NA
  means <- apply(tempEst[,29:33], 1, mean)
  maxs <- apply(tempEst[,29:33], 1, max)
  
  Summary.out <- matrix(nrow=length(means), ncol=3)
  Summary.out[1:length(means), 1] <- unlist(means)
  Summary.out[1:length(means), 2] <- unlist(maxs)
  Summary.out[1:length(means), 3] <- tempEst$RCAID
  
  Summary.out <- as.data.frame(Summary.out)
  colnames(Summary.out) <- c("Mean", "MaxMean", "RCAID")
  sites <- ptIDs[,c("COMID","REACHCODE","Site_ID",paste0(basin, "_rca_id"))]
  ind <- apply(sites, 1, function(x) !any(is.na(x)))
  sites <- sites[ind,]
  
  merged <- merge(sites, Summary.out, by.x=paste0(basin, "_rca_id"), by.y="RCAID", all.x=FALSE, all.y=FALSE)
  
  colnames(merged)[1] <- "RCAID"
  
  merge.mikes <- merge(mikesPts, merged, by.x="SiteName", by.y="Site_ID", all.x=FALSE, all.y=FALSE)
  
  merge.mikes$year <- yearPath
  
  #mikes.out <- merge.mikes  # only once
  
  mikes.out <- rbind(mikes.out, merge.mikes)

  write.csv(mikes.out2, file="Mikes_pts_mn_mxmn.csv", row.names=FALSE)
  
  
  
  
  