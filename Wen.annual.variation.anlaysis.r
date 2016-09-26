############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for annual 8-day means from MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# 
# Edited Aug 2014 to add the PRESS statistic

# Edited Jan 2016 to update the gap-filling interpolation functions
# Edited April 2016 to add model quality testing and pull hard-coded file paths out of the code.

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
library(TeachingDemos)
library(calibrate)
library(rgdal)
library(raster)
library(rNOMADS)
library(gdalUtils)


  basin <- "Wen"
  midBasin <- "Wenatchee"
  longBasin <- "Wenatchee"
  dataType <- "Temp"
  yrPath <- "00"
  monthPath <- "01"
  yearPath <- "2000"
  dataPath <- "D:/OneDrive/work/GIS/NARR/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
    
           
#################################################################
# Data investigation 

# ################################################################
  
 
# ###############################
# PRISM raster processing
# ##################################

  
  
  basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
  pts <- data.frame(rPoints[,1], rPoints[,2])
##################
  # precip
# #################
  

  
  for (j in 1981:2015)
    {  
        year <- j
        yr <- substr(j,3,4)
        setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
        varName <- paste0("Cppt", yr)
        ##### create a list of only the grid files in a directory
        
        allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
        xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
        fileList <- allFiles[!allFiles %in% xmlFiles]
        
        ###### read in one grid to get the structure
        
        rRaster <- raster(fileList[1])
        
        ##### clip the raster to a reasonable PNW extent
        
        ext <- extent(-125, -107, 40, 50)
        rExt <- crop(rRaster, ext)
      
        ##### convert the raster to points and build the data structure for the loop
        
        rPoints <- rasterToPoints(rExt)
        
        pts <- data.frame(rPoints[,1], rPoints[,2])
        data.out <- data.frame(pts, extract(rRaster, pts))
        colnames(data.out) <- c("x", "y", paste0(varName,"001"))
        #rPoints1 <- SpatialPoints(coordinates(rPoints), proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
        
        ##### loop thru all the files in the list and add them iteratively
        
        for (i in 2:365)
          {
          
          r2 <- raster(fileList[i])
          extData <- extract(r2, pts)
          sumData <- extData + data.out[,i+1]
          data.out <- cbind(data.out, sumData)
          namer <- sprintf('%03d', i)
          colnames(data.out)[i+2] <- paste0(varName, namer)
          }
        
        data.out$PtID <- 1:89175
        write.dbf(data.out, file= paste0("SumPpt", year, ".dbf"))
        
        
    }

  
  

######################
# if going back later to process
# #####################
  
  for (j in 1981:2015)
    {
      year <- j
      yr <- substr(j,3,4)
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
      data.in <- read.dbf(paste0("SumPpt", year, ".dbf"))
      data.out <- data.in[data.in$PtID %in% basinPts$PtID,]
      assign(paste0("Ppt_Wen_", year), colMeans(data.out)) 
    }
    
  
  list <- mget(paste0("Ppt_Wen_", 1981:2010))
  all <- do.call(rbind, list)
  
  y30_Ppt_mn <- colMeans(all[,3:367])
  
  
########################
# mean first, then accumulate
# #######################
# ##### 
#   pts <- data.frame(basinPts[,1], basinPts[,2])
#   rPoints1 <- SpatialPoints(coordinates(pts), proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
#   
#   
#   data.out <- data.frame(pts, extract(rRaster, pts))
#   rRaster <- raster(fileList[1])
#   rStack <- extract(r2, rPoints1)
#   
#   for (i in 2:365)
#   {
#     
#     r2 <- raster(fileList[i])
#     extData <- extract(r2, pts)
#     rStack <- cbind(rStack, extData)
#   }
#   
#   basinMn <- colMeans(rStack)
#   basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
#   
#   data.out15_Wen_ppt <- data.out[data.out$PtID %in% basinPts$PtID,]
#     
# #####
# # Subsetting and summarizing by basin ppt
# #####
#   
#   basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
#   
#   data.out15_Wen_ppt <- data.out[data.out$PtID %in% basinPts$PtID,]
#   Wen_15_mn_ppt <- colMeans(data.out15_Wen_ppt)
# ####   
  plot(Ppt_Wen_2000[3:367], pch=16, col="cornflowerblue", cex=.7, main = "Mean cumulative precip, Wenatchee basin, 2000-2015", xlab="Julian Day", ylab="Precip", ylim=c(0,1600))
  points(Ppt_Wen_2001[3:367], pch=16, col="cyan4", cex=.7)
  points(Ppt_Wen_2002[3:367], pch=16, col="khaki", cex=.7)
  points(Ppt_Wen_2003[3:367], pch=16, col="deeppink4", cex=.7)
  points(Ppt_Wen_2004[3:367], pch=16, col="deeppink", cex=.7)
  points(Ppt_Wen_2009[3:367], pch=16, col="darkgrey", cex=.7)
  points(Ppt_Wen_2008[3:367], pch=16, col="lightblue", cex=.7)
  points(Ppt_Wen_2007[3:367], pch=16, col="red", cex=.7)
  points(Ppt_Wen_2011[3:367], pch=16, col="palevioletred", cex=.7)
  points(Ppt_Wen_2006[3:367], pch=16, col="darkgoldenrod1", cex=.7)
  points(Ppt_Wen_2005[3:367], pch=16, col="blueviolet", cex=.7)
  points(Ppt_Wen_2010[3:367], pch=16, col="lightgrey", cex=.7)
  points(Ppt_Wen_2012[3:367], pch=16, col="darkolivegreen4", cex=.7)
  points(Ppt_Wen_2013[3:367], pch=16, col="darkturquoise", cex=.7)
  points(Ppt_Wen_2014[3:367], pch=16, col="chartreuse2", cex=.7)
  points(Ppt_Wen_2015[3:367], pch=15, col="blue2", cex=.7)

  points(y30_Ppt_mn[1:365], pch=16, col="black", cex=.8)
  
   legend(x=grconvertX(c(1.0,1.4), from='npc'),
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)

   y30_precip <- c(mean(Ppt_Wen_1981), mean(Ppt_Wen_1982), mean(Ppt_Wen_1983), mean(Ppt_Wen_1984), mean(Ppt_Wen_1985), mean(Ppt_Wen_1986), mean(Ppt_Wen_1987), mean(Ppt_Wen_1988), mean(Ppt_Wen_1989), mean(Ppt_Wen_1990), mean(Ppt_Wen_1991), mean(Ppt_Wen_1992), mean(Ppt_Wen_1993), mean(Ppt_Wen_1994), mean(Ppt_Wen_1995), mean(Ppt_Wen_1996), mean(Ppt_Wen_1997), mean(Ppt_Wen_1998), mean(Ppt_Wen_1999), mean(Ppt_Wen_2000), mean(Ppt_Wen_2001), mean(Ppt_Wen_2002), mean(Ppt_Wen_2003), mean(Ppt_Wen_2004), mean(Ppt_Wen_2005), mean(Ppt_Wen_2006), mean(Ppt_Wen_2007), mean(Ppt_Wen_2008), mean(Ppt_Wen_2009), mean(Ppt_Wen_2010), mean(Ppt_Wen_2011), mean(Ppt_Wen_2012), mean(Ppt_Wen_2013), mean(Ppt_Wen_2014), mean(Ppt_Wen_2015))
   ####
##################
# Temperature
# #################
  ext <- extent(-125, -107, 40, 50)
  mnNames <- NULL
  
    for (j in 1981:2015)
      {
        year <- j
        yr <- substr(j,3,4)
        setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
        varName <- paste0("Ctmp", yr)
        
        
        ##### create a list of only the grid files in a directory
        
        allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
        xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
        fileList <- allFiles[!allFiles %in% xmlFiles]
        
        ###### read in one grid to get the structure
        
        rRaster <- raster(fileList[1])
        
        ##### clip the raster to a reasonable PNW extent
        
        rExt <- crop(rRaster, ext)
        
        ##### convert the raster to points and build the data structure for the loop
        
        rPoints <- rasterToPoints(rExt)
        
        pts <- data.frame(rPoints[,1], rPoints[,2])
        data.out <- data.frame(pts, extract(rRaster, pts))
        data.out[data.out[,3] < 0, 3] <-  0
        colnames(data.out) <- c("x", "y", paste0(varName,"001"))
        #rPoints1 <- SpatialPoints(coordinates(rPoints), proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
        
        ##### loop thru all the files in the list and add them iteratively
        
        for (i in 2:365)
          {
            
            r2 <- raster(fileList[i])
            extData <- extract(r2, pts)
            extData[extData < 0] <- 0
            sumData <- extData + data.out[,i+1]
            data.out <- cbind(data.out, sumData)
            namer <- sprintf('%03d', i)
            colnames(data.out)[i+2] <- paste0(varName, namer)
          }
        
        data.out$PtID <- 1:89175
        
        write.dbf(data.out, file= paste0("SumTmp", year, ".dbf"))
        
      }
    
############# that one time I needed to reread everything ##########
  
  for (j in 1983:2014)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
      test.in <- read.dbf(paste0("SumTmp", year, ".dbf"))
      data.out_Wen <- test.in[test.in$PtID %in% basinPts$PtID,]
      assign(paste0("Tmp_Wen_", year), colMeans(data.out_Wen))
    }
#################### This part just creates a list of all the vector names and row binds them ##############
  
  list <- mget(paste0("Tmp_Wen_", 1981:2015))
  all <- do.call(rbind, list)
  
  
  
######################
# if going back later to process
# #####################
  
  for (j in 1981:2015)
    {
      year <- j
      yr <- substr(j,3,4)
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
      data.in <- read.dbf(paste0("SumTmp", year, ".dbf"))
      data.out <- data.in[data.in$PtID %in% basinPts$PtID,]
      assign(paste0("Tmp_Wen_", year), colMeans(data.in)) 
    }
  
########################
  
  y30_Tmp_mn <- colMeans(all[,3:367])

  
  ##### Subsetting and summarizing by basin temp

  
  
  plot(1:365, Tmp_Wen_2000[3:367], pch=16, col="cornflowerblue", cex=.7, main = "Mean cumulative temp, Wenatchee basin, 2000-2015", xlab="Julian Day", ylab="Temp", ylim=c(0,3600))
  points(Tmp_Wen_2001[3:367], pch=16, col="cyan4", cex=.7)
  points(Tmp_Wen_2002[3:367], pch=16, col="khaki", cex=.7)
  points(Tmp_Wen_2003[3:367], pch=16, col="deeppink4", cex=.7)
  points(Tmp_Wen_2004[3:367], pch=16, col="deeppink", cex=.7)
  points(Tmp_Wen_2005[3:367], pch=16, col="blueviolet", cex=.7)
  points(Tmp_Wen_2006[3:367], pch=16, col="darkgoldenrod1", cex=.7)
  points(Tmp_Wen_2007[3:367], pch=16, col="red", cex=.7)
  points(Tmp_Wen_2008[3:367], pch=16, col="lightblue", cex=.7)
  points(Tmp_Wen_2009[3:367], pch=16, col="darkgrey", cex=.7)
  points(Tmp_Wen_2010[3:367], pch=16, col="lightgrey", cex=.7)
  points(Tmp_Wen_2011[3:367], pch=16, col="palevioletred", cex=.7)
  points(Tmp_Wen_2012[3:367], pch=16, col="darkolivegreen4", cex=.7)
  points(Tmp_Wen_2013[3:367], pch=16, col="darkturquoise", cex=.7)
  points(Tmp_Wen_2014[3:367], pch=16, col="chartreuse2", cex=.7)
  points(Tmp_Wen_2015[3:367], pch=15, col="blue2", cex=.7)
  
  points(y30_Tmp_mn[1:365], pch=16, col="black", cex=.8)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)


  ##### vector of annual mean temps
  y30_temp <- c(mean(Tmp_Wen_1981), mean(Tmp_Wen_1982), mean(Tmp_Wen_1983), mean(Tmp_Wen_1984), mean(Tmp_Wen_1985), mean(Tmp_Wen_1986), mean(Tmp_Wen_1987), mean(Tmp_Wen_1988), mean(Tmp_Wen_1989), mean(Tmp_Wen_1990), mean(Tmp_Wen_1991), mean(Tmp_Wen_1992), mean(Tmp_Wen_1993), mean(Tmp_Wen_1994), mean(Tmp_Wen_1995), mean(Tmp_Wen_1996), mean(Tmp_Wen_1997), mean(Tmp_Wen_1998), mean(Tmp_Wen_1999), mean(Tmp_Wen_2000), mean(Tmp_Wen_2001), mean(Tmp_Wen_2002), mean(Tmp_Wen_2003), mean(Tmp_Wen_2004), mean(Tmp_Wen_2005), mean(Tmp_Wen_2006), mean(Tmp_Wen_2007), mean(Tmp_Wen_2008), mean(Tmp_Wen_2009), mean(Tmp_Wen_2010), mean(Tmp_Wen_2011), mean(Tmp_Wen_2012), mean(Tmp_Wen_2013), mean(Tmp_Wen_2014), mean(Tmp_Wen_2015))
  #
  Ann_mn_temp <- c(mean(Wen_00_mn), mean(Wen_01_mn), mean(Wen_02_mn), mean(Wen_03_mn), mean(Wen_04_mn), mean(Wen_05_mn), mean(Wen_06_mn), mean(Wen_07_mn), mean(Wen_08_mn), mean(Wen_09_mn), mean(Wen_10_mn), mean(Wen_11_mn), mean(Wen_12_mn), mean(Wen_13_mn), mean(Wen_14_mn), mean(Wen_15_mn))

  col.rainbow <- rainbow(15)
  
  plot(y30_temp[1:30], y30_precip[1:30], pch=21, col="black", bg=col.rainbow, main="Mean sum temp/precip (calendar year)", xlab = "Mean temp", ylab = "Mean precip", cex=1.5, xlim=c(900,1600),)
  points(y30_temp[31:35], y30_precip[31:35], pch=22, col="black", bg=col.rainbow, cex=1.5)
  
  y30_mn_temp <- mean(y30_temp)
  y30_mn_precip <- mean(y30_precip)
  points(y30_mn_temp, y30_mn_precip, pch=8, cex=1.4)
  
  plot(y30_temp[20:35], y30_precip[20:35], pch=c(21, 22, 24), col="black", bg=col.rainbow, main="Mean sum temp/precip (calendar year)", xlab = "Mean temp", ylab = "Mean precip", cex=1.5, xlim=c(900,1600),)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=c(16,15,17), bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=col.rainbow, cex=.65, xpd=NA)
  
#########
# calendar year
#########
  Ann_mn_ppt_cy <- c(mean(Wen_01_mn_ppt[3:367]), mean(Wen_02_mn_ppt[3:367]), mean(Wen_03_mn_ppt[3:367]), mean(Wen_04_mn_ppt[3:367]), mean(Wen_05_mn_ppt[3:367]), mean(Wen_06_mn_ppt[3:367]), mean(Wen_07_mn_ppt[3:367]), mean(Wen_08_mn_ppt[3:367]), mean(Wen_09_mn_ppt[3:367]), mean(Wen_10_mn_ppt[3:367]), mean(Wen_11_mn_ppt[3:367]), mean(Wen_12_mn_ppt[3:367]), mean(Wen_13_mn_ppt[3:367]), mean(Wen_14_mn_ppt[3:367]), mean(Wen_15_mn_ppt[3:367]))
  
  Ann_mn_tmp_cy<- c(mean(Wen_01_mn[3:367]), mean(Wen_02_mn[3:367]), mean(Wen_03_mn[3:367]), mean(Wen_04_mn[3:367]), mean(Wen_05_mn[3:367]), mean(Wen_06_mn[3:367]), mean(Wen_07_mn[3:367]), mean(Wen_08_mn[3:367]), mean(Wen_09_mn[3:367]), mean(Wen_10_mn[3:367]), mean(Wen_11_mn[3:367]), mean(Wen_12_mn[3:367]), mean(Wen_13_mn[3:367]), mean(Wen_14_mn[3:367]), mean(Wen_15_mn[3:367]))
  
  plot(Ann_mn_tmp_cy, Ann_mn_ppt_cy, pch=c(21,22,24)[, col="black", bg=col.rainbow, main="Mean sum temp/precip (calendar year)", xlab = "Mean temp", ylab = "Mean precip", cex=1.5,)
  
  # legend(x=grconvertX(c(1.0,1.4), from='npc'), 
  #        y=grconvertY(c(0.5, 0.9), from='npc'), pch=c(16,15,17), col=col.rainbow, bty="n", legend = c("01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), cex=.75, xpd=NA)

  legend("topright", pch=c(15,16,17), col=col.rainbow, bty="n", legend = c("01", "02", "03", "04", "05", "06", "07", "08","09","10", "11", "12", "13", "14", "15"), cex=.7, pt.cex=.95)
  textxy(clust_out$temp, clust_out$precip, clust_out$group, cex=.8, offset = 1.0, pos=2)
#########
# little chunk of code to add a climate year categorical variable.  This is based on hand-waving eye-ballin'. Change and iterate as needed
##########
  
  for(i in 1:15) 
  {
    if((coeffs.out.all$Ann_ppt_cy[i] >650) & (coeffs.out.all$Ann_tmp_cy[i] >1250)) {
      coeffs.out.all$ClimateYr[i] <- "HW"}
  }
  
######
# more fomal clustering analysis:
# #######
  
  y30mn <- c(y30_mn_temp, y30_mn_precip)
  mydata <- cbind(y30_temp[21:35], y30_precip[21:35])
  
  mydata<- rbind(mydata, y30mn)
  d <- dist(mydata, method = "euclidean")
  fit <- hclust(d, method="ward.D2")
  plot(fit)
  groups <- cutree(fit, k=4)
  rect.hclust(fit, k=4, border="red")
  clust_out <- mydata
  clust_out <- cbind(clust_out, groups)
  colnames(clust_out) <- c("temp", "precip", "group")
                           
  plot(clust_out[,1], clust_out[,2], pch=c(21, 22, 24), col="black", bg=col.rainbow, main="Mean sum temp/precip (calendar year)", xlab = "Mean temp", ylab = "Mean precip", cex=1.5, xlim=c(900,1600),)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=c(16,15,17), bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=col.rainbow, cex=.65, xpd=NA)
  
  clust_out <- as.data.frame(clust_out)
  colnames(clust_out) <- c("temp", "precip", "dist", "group")
  textxy(clust_out$temp, clust_out$precip, clust_out$group, cex=.8, offset = 1.0, pos=2)
  textxy(clust_out$temp, clust_out$precip, rownames(clust_out), cex=.8, offset = 1.0, pos=1)
  
  
  dist <- apply(mydata, 1, function(x) sqrt((x[1]-1172.935)^2 + (x[2]-669.6563)^2)) 
  mydata <- as.data.frame(mydata)
  mydata$year <- rownames(mydata)
  mydata$year[1:15] <- 1:15
  mydata$dist <- dist
####################
  
  plot(Ann_mn_tmp_cy, Ann_mn_ppt_cy, pch=c(21,22,24), col="black", bg=col.rainbow, main="Mean sum temp/precip (calendar year)", xlab = "Mean temp", ylab = "Mean precip", cex=1.5,)
  
  legend("topright", pch=c(15,16,17), col=col.rainbow, bty="n", legend = c("01", "02", "03", "04", "05", "06", "07", "08","09","10", "11", "12", "13", "14", "15"), cex=.7, pt.cex=.95)
  textxy(clust_out$Ann_tmp_cy, clust_out$Ann_ppt_cy, clust_out$ward_dist, cex=.8, offset = 1.0, pos=2)

########################################
# parsing basin-wide data into HUCs
# ##################################
  rm(list=ls(pattern="Ppt_JD_")) #first, get rid of a bunch of memory hogging stuff
  
  pts <- basinPts$PtID
  
  
  for (j in 1981:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
      data.in <- read.dbf(paste0("SumPpt", year, ".dbf"))
      data.out <- data.in[data.in$PtID %in% pts,]
      write.dbf(data.out, paste0("SumPpt", basin, year, ".dbf"))
      assign(paste0("Ppt_", basin, "_", year), colMeans(data.in)
    }
   
  
  
  lister <- levels(basinPts$HUC_10)
  
  i <- lister[7]
  HUCpts <- basinPts[basinPts$HUC_10 == i, "PtID"]
  
  for (j in 1981:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
      data.in <- read.dbf(paste0("SumPpt", basin, year, ".dbf"))
      data.out <- data.in[data.in$PtID %in% HUCpts,]
      assign(paste0("Ppt_", basin, "_HUC_", 7, "_", year), colMeans(data.out))
    }
  
  
  listed <- mget(paste0("Ppt_", basin, "_HUC_", 5, "_", 1981:2010))
  all <- do.call(rbind, listed)
  assign(paste0("y30_Ppt_HUC", 5), colMeans(all))
  
  y30_Tmp_HUC5_mn <- mean(y30_Tmp_HUC5[3:367])
  y30_Ppt_HUC5_mn <- mean(y30_Ppt_HUC5[3:367])
  
  Ann_mn_ppt_HUC1 <- c(mean(Ppt_Wen_HUC_1_2000[3:367]), mean(Ppt_Wen_HUC_1_2001[3:367]), mean(Ppt_Wen_HUC_1_2002[3:367]), mean(Ppt_Wen_HUC_1_2003[3:367]), mean(Ppt_Wen_HUC_1_2004[3:367]), mean(Ppt_Wen_HUC_1_2005[3:367]), mean(Ppt_Wen_HUC_1_2006[3:367]), mean(Ppt_Wen_HUC_1_2007[3:367]), mean(Ppt_Wen_HUC_1_2008[3:367]), mean(Ppt_Wen_HUC_1_2009[3:367]), mean(Ppt_Wen_HUC_1_2010[3:367]), mean(Ppt_Wen_HUC_1_2011[3:367]), mean(Ppt_Wen_HUC_1_2012[3:367]), mean(Ppt_Wen_HUC_1_2013[3:367]), mean(Ppt_Wen_HUC_1_2014[3:367]), mean(Ppt_Wen_HUC_1_2015[3:367]))
  
  Ann_mn_tmp_HUC1<- c(mean(Tmp_Wen_HUC_1_2000[3:367]), mean(Tmp_Wen_HUC_1_2001[3:367]), mean(Tmp_Wen_HUC_1_2002[3:367]), mean(Tmp_Wen_HUC_1_2003[3:367]), mean(Tmp_Wen_HUC_1_2004[3:367]), mean(Tmp_Wen_HUC_1_2005[3:367]), mean(Tmp_Wen_HUC_1_2006[3:367]), mean(Tmp_Wen_HUC_1_2007[3:367]), mean(Tmp_Wen_HUC_1_2008[3:367]), mean(Tmp_Wen_HUC_1_2009[3:367]), mean(Tmp_Wen_HUC_1_2010[3:367]), mean(Tmp_Wen_HUC_1_2011[3:367]), mean(Tmp_Wen_HUC_1_2012[3:367]), mean(Tmp_Wen_HUC_1_2013[3:367]), mean(Tmp_Wen_HUC_1_2014[3:367]), mean(Tmp_Wen_HUC_1_2015[3:367]))
  
  y30_temp <- c(mean(Tmp_Wen_1981), mean(Tmp_Wen_1982), mean(Tmp_Wen_1983), mean(Tmp_Wen_1984), mean(Tmp_Wen_1985), mean(Tmp_Wen_1986), mean(Tmp_Wen_1987), mean(Tmp_Wen_1988), mean(Tmp_Wen_1989), mean(Tmp_Wen_1990), mean(Tmp_Wen_1991), mean(Tmp_Wen_1992), mean(Tmp_Wen_1993), mean(Tmp_Wen_1994), mean(Tmp_Wen_1995), mean(Tmp_Wen_1996), mean(Tmp_Wen_1997), mean(Tmp_Wen_1998), mean(Tmp_Wen_1999), mean(Tmp_Wen_2000), mean(Tmp_Wen_2001), mean(Tmp_Wen_2002), mean(Tmp_Wen_2003), mean(Tmp_Wen_2004), mean(Tmp_Wen_2005), mean(Tmp_Wen_2006), mean(Tmp_Wen_2007), mean(Tmp_Wen_2008), mean(Tmp_Wen_2009), mean(Tmp_Wen_2010), mean(Tmp_Wen_2011), mean(Tmp_Wen_2012), mean(Tmp_Wen_2013), mean(Tmp_Wen_2014), mean(Tmp_Wen_2015))
  
  y01_15_temp <- c(mean(Tmp_Wen_2001), mean(Tmp_Wen_2002), mean(Tmp_Wen_2003), mean(Tmp_Wen_2004), mean(Tmp_Wen_2005), mean(Tmp_Wen_2006), mean(Tmp_Wen_2007), mean(Tmp_Wen_2008), mean(Tmp_Wen_2009), mean(Tmp_Wen_2010), mean(Tmp_Wen_2011), mean(Tmp_Wen_2012), mean(Tmp_Wen_2013), mean(Tmp_Wen_2014), mean(Tmp_Wen_2015))
  
  y30_ppt <- c(mean(Ppt_Wen_1981), mean(Ppt_Wen_1982), mean(Ppt_Wen_1983), mean(Ppt_Wen_1984), mean(Ppt_Wen_1985), mean(Ppt_Wen_1986), mean(Ppt_Wen_1987), mean(Ppt_Wen_1988), mean(Ppt_Wen_1989), mean(Ppt_Wen_1990), mean(Ppt_Wen_1991), mean(Ppt_Wen_1992), mean(Ppt_Wen_1993), mean(Ppt_Wen_1994), mean(Ppt_Wen_1995), mean(Ppt_Wen_1996), mean(Ppt_Wen_1997), mean(Ppt_Wen_1998), mean(Ppt_Wen_1999), mean(Ppt_Wen_2000), mean(Ppt_Wen_2001), mean(Ppt_Wen_2002), mean(Ppt_Wen_2003), mean(Ppt_Wen_2004), mean(Ppt_Wen_2005), mean(Ppt_Wen_2006), mean(Ppt_Wen_2007), mean(Ppt_Wen_2008), mean(Ppt_Wen_2009), mean(Ppt_Wen_2010), mean(Ppt_Wen_2011), mean(Ppt_Wen_2012), mean(Ppt_Wen_2013), mean(Ppt_Wen_2014), mean(Ppt_Wen_2015))
  
  y01_15_ppt <- c(mean(Ppt_Wen_2001), mean(Ppt_Wen_2002), mean(Ppt_Wen_2003), mean(Ppt_Wen_2004), mean(Ppt_Wen_2005), mean(Ppt_Wen_2006), mean(Ppt_Wen_2007), mean(Ppt_Wen_2008), mean(Ppt_Wen_2009), mean(Ppt_Wen_2010), mean(Ppt_Wen_2011), mean(Ppt_Wen_2012), mean(Ppt_Wen_2013), mean(Ppt_Wen_2014), mean(Ppt_Wen_2015))
  
  
  Ann_mn_HUC_7 <- data.frame(Ppt=Ann_mn_ppt_HUC_7, Tmp=Ann_mn_tmp_HUC_7, year=sprintf('%02d', 0:15), HUC=7)
  
  list <- mget(paste0("Ppt_", basin, "_", 1981:2010))
  all <- do.call(rbind, list)
  y30_precip <- colMeans(all[,3:367])
  y30_precip_mn <- mean(y30_precip)
  
  plot(Ann_mn_tmp_HUC1, Ann_mn_ppt_HUC1, pch=21, col="black", bg="lightblue", main="Mean sum temp/precip Wen (by HUC) 2000-2015", xlab = "Mean temp", ylab = "Mean precip", cex=1.5, xlim=c(600,1900), ylim=c(100,1300),)
  points(Ann_mn_tmp_HUC2, Ann_mn_ppt_HUC2, pch=22, col="black", bg="orange", cex=1.5)
  points(Ann_mn_tmp_HUC_3, Ann_mn_ppt_HUC_3, pch=24, col="black", bg="lightgreen", cex=1.5)
  points(Ann_mn_tmp_HUC_4, Ann_mn_ppt_HUC_4, pch=23, col="black", bg="yellow", cex=1.5)
  points(Ann_mn_tmp_HUC_5, Ann_mn_ppt_HUC_5, pch=1, col="blue", bg="lightblue", cex=1.5)
  points(Ann_mn_tmp_HUC_6, Ann_mn_ppt_HUC_6, pch=2, col="green", bg="lightgreen", cex=1.5)
  points(Ann_mn_tmp_HUC_7, Ann_mn_ppt_HUC_7, pch=7, col="red", bg="orange", cex=1.5)
  points(y30_Tmp_HUC1_mn, y30_Ppt_HUC1_mn, pch="#", col="lightblue", cex=1.5)
  points(y30_Tmp_HUC2_mn, y30_Ppt_HUC2_mn, pch="#", col="orange", cex=1.5)
  points(y30_Tmp_HUC3_mn, y30_Ppt_HUC3_mn, pch="#", col="lightgreen", cex=1.5)
  points(y30_Tmp_HUC4_mn, y30_Ppt_HUC4_mn, pch="#", col="yellow", cex=1.5)
  points(y30_Tmp_HUC5_mn, y30_Ppt_HUC5_mn, pch="#", col="blue", cex=1.5)
  points(y30_Tmp_HUC6_mn, y30_Ppt_HUC6_mn, pch="#", col="green", cex=1.5)
  points(y30_Tmp_HUC7_mn, y30_Ppt_HUC7_mn, pch="#", col="red", cex=1.5)
  points(y30_Tmp_mn, y30_precip_mn, pch=8, col="black", cex=2.0)
  
  legend("topright", pch=c(21,22,24,23,1,2,7,8), col=c("black", "black", "lightgreen", "yellow", "blue", "green", "red", "black"), bg=c("lightblue", "orange", "lightgreen", "yellow", "blue", "green", "red", "black"),bty="n", legend = c("SFJD", "NFJD", "MFJD", "LFJD", "30yJD"), cex=.7, pt.cex=.95)
  textxy(Ann_mn_HUC_1$Tmp, Ann_mn_HUC_1$Ppt, Ann_mn_HUC_1$year, cex=.8, offset = 1.0, pos=2)
  
  listed <- mget(paste0("Ppt_", basin, "_HUC_", 5, "_", 1981:2015))
  all <- do.call(rbind, listed)
  write.csv(all, paste0("Ppt_", basin, "_HUC_1981_2015.csv"))
###########################################
# parsing annual data into water-year sets
##########
  
  ext <- extent(-125, -107, 40, 50)
  basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
  pts <- data.frame(rPoints[,1], rPoints[,2])
  dirPath <- "D:/OneDrive/work/GIS/PRISM/"
  
  for (k in 10:15)
    {
      year1 <- paste0("200", k)
      yr1 <- paste0("0", k)
      year2 <- paste0("200", k+1)
      yr2 <- paste0("0", k+1)
      
        setwd(paste0(dirPath, year1))
        varName <- paste0("Cppt", yr1)
        
        ##### create lists of only the grid files in 2 annual directoies
        
        allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
        xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
        fileList1 <- allFiles[!allFiles %in% xmlFiles]
        rRaster1 <- raster(fileList1[244])
        
        setwd(paste0(dirPath, year2))
        allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
        xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
        fileList2 <- allFiles[!allFiles %in% xmlFiles]
        
        
        ##### clip the rasters to a reasonable PNW extent
        
        rExt <- crop(rRaster1, ext)
        
        ##### convert the rasters to points and build the data structure for the loop
        
        rPoints <- rasterToPoints(rExt)
        
        data.out <- data.frame(pts, extract(rRaster, pts))
        data.out[data.out[,3] < 0, 3] <-  0
        colnames(data.out) <- c("x", "y", paste0(varName,"244"))
        
        ##### loop thru all the files in the lists and add them iteratively
        
        setwd(paste0(dirPath, year1))
        j <- 2
        
        for (i in 245:365)
          {
            r2 <- raster(fileList1[i])
            extData <- extract(r2, pts)
            extData[extData < 0] <- 0
            sumData <- extData + data.out[,j+1]
            data.out <- cbind(data.out, sumData)
            namer <- sprintf('%03d', i)
            colnames(data.out)[j+2] <- paste0(varName, namer)
            j <- j+1
          }
        
        setwd(paste0(dirPath, year2))
        varName <- paste0("Cppt", yr2)
        
        for (i in 1:243)
          {
            
            r2 <- raster(fileList2[i])
            extData <- extract(r2, pts)
            extData[extData < 0] <- 0
            sumData <- extData + data.out[,j+1]
            data.out <- cbind(data.out, sumData)
            namer <- sprintf('%03d', i)
            colnames(data.out)[j+2] <- paste0(varName, namer)
            j <- j+1
          }
        
        data.out$PtID <- 1:89175
        data.out[3:367,] <- round(data.out[3:367,], 2) 
        setwd("D:/OneDrive/work/GIS/PRISM/water_year_output/")
        write.dbf(data.out, file= paste0("SumPpt", yr1, yr2, ".dbf"))
        
        data.out_Wen <- data.out[data.out$PtID %in% basinPts$PtID,]
        #assign(paste0("Wen_", yr1, yr2, "_mn"), colMeans(data.out_Wen))
        assign(paste0("Wen_", yr1, yr2, "_mn_ppt"), colMeans(data.out_Wen))
        
    }
  
  
#####
# Subsetting and summarizing by basin
#####
  
  
  plot(Wen_0708_mn, pch=16, col="lightgrey", cex=.7, xaxt="n", main = "Mean cumulative temp by water-year, Wen 2000-2015", xlab="Sept-Aug", ylab="Temp", ylim=c(0,3400))
  points(Wen_0809_mn, pch=16, col="black", cex=.7)
  points(Wen_0910_mn, pch=16, col="blue2", cex=.7)
  points(Wen_1011_mn, pch=16, col="red", cex=.7)
  points(Wen_1112_mn, pch=16, col="palevioletred", cex=.7)
  points(Wen_1213_mn, pch=16, col="darkolivegreen4", cex=.7)
  points(Wen_1314_mn, pch=16, col="darkturquoise", cex=.7)
  points(Wen_1415_mn, pch=16, col="chartreuse2", cex=.7)
  points(Wen_0001_mn, pch=16, col="cornflowerblue", cex=.7)
  points(Wen_0102_mn, pch=16, col="cyan4", cex=.7)
  points(Wen_0203_mn, pch=16, col="khaki", cex=.7)
  points(Wen_0304_mn, pch=16, col="deeppink4", cex=.7)
  points(Wen_0405_mn, pch=16, col="deeppink", cex=.7)
  points(Wen_0506_mn, pch=16, col="blueviolet ", cex=.7)
  points(Wen_0607_mn, pch=16, col="darkgoldenrod1  ", cex=.7)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708","0809","0910", "1011", "1112", "1213", "1314", "1415"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "goldenrod1", "lightgrey", "black", "blue2", "red", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2"), cex=.6, xpd=NA)

  temp_pal <- c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "goldenrod1", "lightgrey", "black", "blue2", "red", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2")
  col.rainbow <- rainbow(15)
  
  Ann_mn_temp_wy <- c(mean(Wen_0001_mn), mean(Wen_0102_mn), mean(Wen_0203_mn), mean(Wen_0304_mn), mean(Wen_0405_mn), mean(Wen_0506_mn), mean(Wen_0607_mn), mean(Wen_0708_mn), mean(Wen_0809_mn), mean(Wen_0910_mn), mean(Wen_1011_mn), mean(Wen_1112_mn), mean(Wen_1213_mn), mean(Wen_1314_mn), mean(Wen_1415_mn))
  Ann_mn_ppt_wy <- c(mean(Wen_0001_mn_ppt), mean(Wen_0102_mn_ppt), mean(Wen_0203_mn_ppt), mean(Wen_0304_mn_ppt), mean(Wen_0405_mn_ppt), mean(Wen_0506_mn_ppt), mean(Wen_0607_mn_ppt), mean(Wen_0708_mn_ppt), mean(Wen_0809_mn_ppt), mean(Wen_0910_mn_ppt), mean(Wen_1011_mn_ppt), mean(Wen_1112_mn_ppt), mean(Wen_1213_mn_ppt), mean(Wen_1314_mn_ppt), mean(Wen_1415_mn_ppt))
  
  plot(Ann_mn_temp_wy, Ann_mn_ppt_wy, pch=c(21,22,24), col="black", bg=col.rainbow, main="Mean sum temp/precip by water year", xlab = "Annual mean temp", ylab = "Annual mean precip", cex=1.5)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=c(21,22,23), col="black", bg=col.rainbow, bty="n", legend = c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708","0809","0910", "1011", "1112", "1213", "1314", "1415"), cex=.7, xpd=NA)
 
  legend("topright", pch=c(16, 15, 17), col=col.rainbow, bty="n", legend = c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708","0809","0910", "1011", "1112", "1213", "1314", "1415"), cex=.7, pt.cex=1.0)
  
  
#####################
# plot temp by precip
####################
  
  plot(Wen_00_mn_ppt[3:367], Wen_00_mn[3:367], pch=16, col="cornflowerblue", cex=.7, main = "Mean cumulative temp & precip, Wen basin, 2000-2015", xlab="Precip", ylab="Temp", ylim=c(-50,1700))
  points(Wen_01_mn_ppt[3:367], pch=16, col="cyan4", cex=.7)
  points(Wen_02_mn_ppt[3:367], pch=16, col="khaki", cex=.7)
  points(Wen_03_mn_ppt[3:367], pch=16, col="deeppink4", cex=.7)
  points(Wen_04_mn_ppt[3:367], pch=16, col="deeppink", cex=.7)
  points(Wen_09_mn_ppt[3:367], pch=16, col="black", cex=.7)
  points(Wen_08_mn_ppt[3:367], pch=16, col="lightblue", cex=.7)
  points(Wen_07_mn_ppt[3:367], pch=16, col="red", cex=.7)
  points(Wen_11_mn_ppt[3:367], pch=16, col="palevioletred", cex=.7)
  points(Wen_06_mn_ppt[3:367], pch=16, col="darkgoldenrod1", cex=.7)
  points(Wen_05_mn_ppt[3:367], pch=16, col="blueviolet", cex=.7)
  points(Wen_10_mn_ppt[3:367], pch=16, col="lightgrey", cex=.7)
  points(Wen_12_mn_ppt[3:367], pch=16, col="darkolivegreen4", cex=.7)
  points(Wen_13_mn_ppt[3:367], pch=16, col="darkturquoise", cex=.7)
  points(Wen_14_mn_ppt[3:367], pch=16, col="chartreuse2", cex=.7)
  points(Wen_15_mn_ppt[3:367], pch=15, col="blue2", cex=.7)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)
 
###########################################
# parsing annual data into spring/fall sets
# precip
##########
 
  ext <- extent(-125, -107, 40, 50)
  basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
  pts <- data.frame(rPoints[,1], rPoints[,2])
  dirPath <- "D:/OneDrive/work/GIS/PRISM/"
  
  for (k in 10:15)
      {
        year <- paste0("20", k)
        yr <- paste0("", k)
        
        setwd(paste0(dirPath, year))
        varName <- paste0("Cppt", yr)
        
        ##### create lists of only the grid files in 2 annual directoies
        
        allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
        xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
        fileList <- allFiles[!allFiles %in% xmlFiles]
        rRaster <- raster(fileList[1])
        
        ##### clip the rasters to a reasonable PNW extent
        
        rExt <- crop(rRaster, ext)
        
        ##### convert the rasters to points and build the data structure for the loop
        
        rPoints <- rasterToPoints(rExt)
        
        data.out <- data.frame(pts, extract(rRaster, pts))
        data.out[data.out[,3] < 0, 3] <-  0
        colnames(data.out) <- c("x", "y", paste0(varName,"001"))
        
        ##### loop thru all the files in the lists and add them iteratively
        
        setwd(paste0(dirPath, year))
        j <- 2
        seas <- "sp"
        
        for (i in 2:181)
        {
          r2 <- raster(fileList[i])
          extData <- extract(r2, pts)
          extData[extData < 0] <- 0
          sumData <- extData + data.out[,j+1]
          data.out <- cbind(data.out, sumData)
          namer <- sprintf('%03d', i)
          colnames(data.out)[j+2] <- paste0(varName, namer)
          j <- j+1
        }
        
        data.out$PtID <- 1:89175
        data.out[3:183,] <- round(data.out[3:183,], 2) 
        setwd("D:/OneDrive/work/GIS/PRISM/water_year_output/halfYears")
        write.dbf(data.out, file= paste0("SumPpt", yr, seas, ".dbf"))
        
        data.out_Wen <- data.out[data.out$PtID %in% basinPts$PtID,]
        
        assign(paste0("Wen_", yr, seas, "_mn_ppt"), colMeans(data.out_Wen))
        
        setwd(paste0(dirPath, year))
        j <- 2
        rRaster <- raster(fileList[182])
        data.out <- data.frame(pts, extract(rRaster, pts))
        data.out[data.out[,3] < 0, 3] <-  0
        colnames(data.out) <- c("x", "y", paste0(varName,"182"))
        
        seas2 <- "fall"
        
        for (i in 182:365)
        {
          
          r2 <- raster(fileList[i])
          extData <- extract(r2, pts)
          extData[extData < 0] <- 0
          sumData <- extData + data.out[,j+1]
          data.out <- cbind(data.out, sumData)
          namer <- sprintf('%03d', i)
          colnames(data.out)[j+2] <- paste0(varName, namer)
          j <- j+1
        }
        
        data.out$PtID <- 1:89175
        data.out[3:187,] <- round(data.out[3:187,], 2) 
        setwd("D:/OneDrive/work/GIS/PRISM/water_year_output/halfYears")
        write.dbf(data.out, file= paste0("SumPpt", yr, seas2, ".dbf"))
        
        data.out_Wen <- data.out[data.out$PtID %in% basinPts$PtID,]
        
        assign(paste0("Wen_", yr, seas2, "_mn_ppt"), colMeans(data.out_Wen))
        
      }
      
  
#####
# Plot & summarizing 
#####
  
  
  plot(Wen_00sp_mn_ppt[3:183], pch=21, col="black", bg="cornflowerblue", cex=.9, main = "Mean cumulative precip sp, Wenatchee basin, 2000-2015", xlab="Julian Day", ylab="Precip", ylim=c(0,1000))
  points(Wen_00fall_mn_ppt[3:187], pch=17, col="cornflowerblue", cex=.9)
  points(Wen_01sp_mn_ppt[3:183], pch=21, col="black", bg="cyan4", cex=.9)
  points(Wen_01fall_mn_ppt[3:187], pch=17, col="cyan4", cex=.9)
  points(Wen_02sp_mn_ppt[3:183], pch=21, col="black", bg="khaki", cex=.9)
  points(Wen_02fall_mn_ppt[3:187], pch=17, col="khaki", cex=.9)
  points(Wen_03sp_mn_ppt[3:183], pch=21, col="black", bg="deeppink4", cex=.9)
  points(Wen_03fall_mn_ppt[3:187], pch=17, col="deeppink4", cex=.9)
  points(Wen_04sp_mn_ppt[3:183], pch=21, col="black", bg="deeppink", cex=.9)
  points(Wen_04fall_mn_ppt[3:187], pch=17, col="deeppink", cex=.9)
  points(Wen_05sp_mn_ppt[3:183], pch=21, col="black", bg="blueviolet", cex=.9)
  points(Wen_05fall_mn_ppt[3:187], pch=17, col="blueviolet", cex=.9)
  points(Wen_06sp_mn_ppt[3:183], pch=21, col="black", bg="darkgoldenrod1", cex=.9)
  points(Wen_06fall_mn_ppt[3:187], pch=17, col="darkgoldenrod1", cex=.9)
  points(Wen_07sp_mn_ppt[3:183], pch=21, col="black", bg="red", cex=.9)
  points(Wen_07fall_mn_ppt[3:187], pch=17, col="red", cex=.9)
  points(Wen_08sp_mn_ppt[3:183], pch=21, col="black", bg="lightblue", cex=.9)
  points(Wen_08fall_mn_ppt[3:187], pch=17, col="lightblue", cex=.9)
  points(Wen_09sp_mn_ppt[3:183], pch=21, col="black", bg="black", cex=.9)
  points(Wen_09fall_mn_ppt[3:187], pch=17, col="black", cex=.9)
  points(Wen_10sp_mn_ppt[3:183], pch=21, col="black", bg="lightgrey", cex=.9)
  points(Wen_10fall_mn_ppt[3:187], pch=17, col="lightgrey", cex=.9)
  points(Wen_11sp_mn_ppt[3:183], pch=21, col="black", bg="palevioletred", cex=.9)
  points(Wen_11fall_mn_ppt[3:187], pch=17, col="palevioletred", cex=.9)
  points(Wen_12sp_mn_ppt[3:183], pch=21, col="black", bg="darkolivegreen4", cex=.9)
  points(Wen_12fall_mn_ppt[3:187], pch=17, col="darkolivegreen4", cex=.9)
  points(Wen_13sp_mn_ppt[3:183], pch=21, col="black", bg="darkturquoise", cex=.9)
  points(Wen_13fall_mn_ppt[3:187], pch=17, col="darkturquoise", cex=.9)
  points(Wen_14sp_mn_ppt[3:183], pch=21, col="black", bg="chartreuse2", cex=.9)
  points(Wen_14fall_mn_ppt[3:187], pch=17, col="chartreuse2", cex=.9)
  points(Wen_15sp_mn_ppt[3:183], pch=21, col="black", bg="blue2", cex=.9)
  points(Wen_15fall_mn_ppt[3:187], pch=17, col="blue2", cex=.9)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)
  
  
  temp_pal <- c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "goldenrod1", "lightgrey", "black", "blue2", "red", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2")
  col.rainbow <- rainbow(15)
  
  Ann_mn_ppt_sp <- c(mean(Wen_01sp_mn_ppt), mean(Wen_02sp_mn_ppt), mean(Wen_03sp_mn_ppt), mean(Wen_04sp_mn_ppt), mean(Wen_05sp_mn_ppt), mean(Wen_06sp_mn_ppt), mean(Wen_07sp_mn_ppt), mean(Wen_08sp_mn_ppt), mean(Wen_09sp_mn_ppt), mean(Wen_10sp_mn_ppt), mean(Wen_11sp_mn_ppt), mean(Wen_12sp_mn_ppt), mean(Wen_13sp_mn_ppt), mean(Wen_14sp_mn_ppt), mean(Wen_15sp_mn_ppt))
  Ann_mn_ppt_fall <- c(mean(Wen_01fall_mn_ppt), mean(Wen_02fall_mn_ppt), mean(Wen_03fall_mn_ppt), mean(Wen_04fall_mn_ppt), mean(Wen_05fall_mn_ppt), mean(Wen_06fall_mn_ppt), mean(Wen_07fall_mn_ppt), mean(Wen_08fall_mn_ppt), mean(Wen_09fall_mn_ppt), mean(Wen_10fall_mn_ppt), mean(Wen_11fall_mn_ppt), mean(Wen_12fall_mn_ppt), mean(Wen_13fall_mn_ppt), mean(Wen_14fall_mn_ppt), mean(Wen_15fall_mn_ppt))
  
  plot(Ann_mn_temp_wy, Ann_mn_ppt_wy, pch=c(21,22,24), col="black", bg=col.rainbow, main="Mean sum temp/precip by water year", xlab = "Annual mean temp", ylab = "Annual mean precip", cex=1.5)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=c(21,22,23), col="black", bg=col.rainbow, bty="n", legend = c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708","0809","0910", "1011", "1112", "1213", "1314", "1415"), cex=.7, xpd=NA)
  
  legend("topright", pch=c(16, 15, 17), col=col.rainbow, bty="n", legend = c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708","0809","0910", "1011", "1112", "1213", "1314", "1415"), cex=.7, pt.cex=1.0)
  
###########################################
# parsing annual data into spring/fall sets
# Temp
##########
  
  ext <- extent(-125, -107, 40, 50)
  basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
  pts <- data.frame(rPoints[,1], rPoints[,2])
  dirPath <- "D:/OneDrive/work/GIS/PRISM/Temp/"
  
  for (k in 10:15)
    {
      year <- paste0("20", k)
      yr <- paste0("", k)
      
      setwd(paste0(dirPath, year))
      varName <- paste0("Ctmp", yr)
      
      ##### create lists of only the grid files in 2 annual directoies
      
      allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
      xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
      fileList <- allFiles[!allFiles %in% xmlFiles]
      rRaster <- raster(fileList[1])
      
      ##### clip the rasters to a reasonable PNW extent
      
      rExt <- crop(rRaster, ext)
      
      ##### convert the rasters to points and build the data structure for the loop
      
      rPoints <- rasterToPoints(rExt)
      
      data.out <- data.frame(pts, extract(rRaster, pts))
      data.out[data.out[,3] < 0, 3] <-  0
      colnames(data.out) <- c("x", "y", paste0(varName,"001"))
      
      ##### loop thru all the files in the lists and add them iteratively
      
      setwd(paste0(dirPath, year))
      j <- 2
      seas <- "sp"
      
      for (i in 2:181)
        {
          r2 <- raster(fileList[i])
          extData <- extract(r2, pts)
          extData[extData < 0] <- 0
          sumData <- extData + data.out[,j+1]
          data.out <- cbind(data.out, sumData)
          namer <- sprintf('%03d', i)
          colnames(data.out)[j+2] <- paste0(varName, namer)
          j <- j+1
        }
      
      data.out$PtID <- 1:89175
      data.out[3:183,] <- round(data.out[3:183,], 2) 
      setwd("D:/OneDrive/work/GIS/PRISM/water_year_output/halfYears")
      write.dbf(data.out, file= paste0("SumTmp", yr, seas, ".dbf"))
      
      data.out_Wen <- data.out[data.out$PtID %in% basinPts$PtID,]
      
      assign(paste0("Wen_", yr, seas, "_mn_tmp"), colMeans(data.out_Wen))
      
      setwd(paste0(dirPath, year))
      j <- 2
      rRaster <- raster(fileList[182])
      data.out <- data.frame(pts, extract(rRaster, pts))
      data.out[data.out[,3] < 0, 3] <-  0
      colnames(data.out) <- c("x", "y", paste0(varName,"182"))
      
      seas2 <- "fall"
      
      for (i in 182:365)
        {
          
          r2 <- raster(fileList[i])
          extData <- extract(r2, pts)
          extData[extData < 0] <- 0
          sumData <- extData + data.out[,j+1]
          data.out <- cbind(data.out, sumData)
          namer <- sprintf('%03d', i)
          colnames(data.out)[j+2] <- paste0(varName, namer)
          j <- j+1
        }
      
      data.out$PtID <- 1:89175
      data.out[3:187,] <- round(data.out[3:187,], 2) 
      setwd("D:/OneDrive/work/GIS/PRISM/water_year_output/halfYears")
      write.dbf(data.out, file= paste0("SumTmp", yr, seas2, ".dbf"))
      
      data.out_Wen <- data.out[data.out$PtID %in% basinPts$PtID,]
      
      assign(paste0("Wen_", yr, seas2, "_mn_tmp"), colMeans(data.out_Wen))
      
    }
  
  
#####
# Plot & summarizing 
#####
  
  
  plot(Wen_00sp_mn_tmp[3:183], pch=21, col="black", bg="cornflowerblue", cex=.9, main = "Mean cumulative temp sp, Wenatchee basin, 2000-2015", xlab="Julian Day", ylab="Temp", ylim=c(0,1500))
  points(Wen_00fall_mn_tmp[3:187], pch=17, col="cornflowerblue", cex=.9)
  points(Wen_01sp_mn_tmp[3:183], pch=21, col="black", bg="cyan4", cex=.9)
  points(Wen_01fall_mn_tmp[3:187], pch=17, col="cyan4", cex=.9)
  points(Wen_02sp_mn_tmp[3:183], pch=21, col="black", bg="khaki", cex=.9)
  points(Wen_02fall_mn_tmp[3:187], pch=17, col="khaki", cex=.9)
  points(Wen_03sp_mn_tmp[3:183], pch=21, col="black", bg="deeppink4", cex=.9)
  points(Wen_03fall_mn_tmp[3:187], pch=17, col="deeppink4", cex=.9)
  points(Wen_04sp_mn_tmp[3:183], pch=21, col="black", bg="deeppink", cex=.9)
  points(Wen_04fall_mn_tmp[3:187], pch=17, col="deeppink", cex=.9)
  points(Wen_05sp_mn_tmp[3:183], pch=21, col="black", bg="blueviolet", cex=.9)
  points(Wen_05fall_mn_tmp[3:187], pch=17, col="blueviolet", cex=.9)
  points(Wen_06sp_mn_tmp[3:183], pch=21, col="black", bg="darkgoldenrod1", cex=.9)
  points(Wen_06fall_mn_tmp[3:187], pch=17, col="darkgoldenrod1", cex=.9)
  points(Wen_07sp_mn_tmp[3:183], pch=21, col="black", bg="red", cex=.9)
  points(Wen_07fall_mn_tmp[3:187], pch=17, col="red", cex=.9)
  points(Wen_08sp_mn_tmp[3:183], pch=21, col="black", bg="lightblue", cex=.9)
  points(Wen_08fall_mn_tmp[3:187], pch=17, col="lightblue", cex=.9)
  points(Wen_09sp_mn_tmp[3:183], pch=21, col="black", bg="black", cex=.9)
  points(Wen_09fall_mn_tmp[3:187], pch=17, col="black", cex=.9)
  points(Wen_10sp_mn_tmp[3:183], pch=21, col="black", bg="lightgrey", cex=.9)
  points(Wen_10fall_mn_tmp[3:187], pch=17, col="lightgrey", cex=.9)
  points(Wen_11sp_mn_tmp[3:183], pch=21, col="black", bg="palevioletred", cex=.9)
  points(Wen_11fall_mn_tmp[3:187], pch=17, col="palevioletred", cex=.9)
  points(Wen_12sp_mn_tmp[3:183], pch=21, col="black", bg="darkolivegreen4", cex=.9)
  points(Wen_12fall_mn_tmp[3:187], pch=17, col="darkolivegreen4", cex=.9)
  points(Wen_13sp_mn_tmp[3:183], pch=21, col="black", bg="darkturquoise", cex=.9)
  points(Wen_13fall_mn_tmp[3:187], pch=17, col="darkturquoise", cex=.9)
  points(Wen_14sp_mn_tmp[3:183], pch=21, col="black", bg="chartreuse2", cex=.9)
  points(Wen_14fall_mn_tmp[3:187], pch=17, col="chartreuse2", cex=.9)
  points(Wen_15sp_mn_tmp[3:183], pch=21, col="black", bg="blue2", cex=.9)
  points(Wen_15fall_mn_tmp[3:187], pch=17, col="blue2", cex=.9)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)
  
  
  temp_pal <- c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "goldenrod1", "lightgrey", "black", "blue2", "red", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2")
  col.rainbow <- rainbow(15)
  
  Ann_mn_ppt_sp <- c(mean(Wen_01sp_mn_ppt[3:183]), mean(Wen_02sp_mn_ppt[3:183]), mean(Wen_03sp_mn_ppt[3:183]), mean(Wen_04sp_mn_ppt[3:183]), mean(Wen_05sp_mn_ppt[3:183]), mean(Wen_06sp_mn_ppt[3:183]), mean(Wen_07sp_mn_ppt[3:183]), mean(Wen_08sp_mn_ppt[3:183]), mean(Wen_09sp_mn_ppt[3:183]), mean(Wen_10sp_mn_ppt[3:183]), mean(Wen_11sp_mn_ppt[3:183]), mean(Wen_12sp_mn_ppt[3:183]), mean(Wen_13sp_mn_ppt[3:183]), mean(Wen_14sp_mn_ppt[3:183]), mean(Wen_15sp_mn_ppt[3:183]))
  Ann_mn_ppt_fall <- c(mean(Wen_01fall_mn_ppt[3:187]), mean(Wen_02fall_mn_ppt[3:187]), mean(Wen_03fall_mn_ppt[3:187]), mean(Wen_04fall_mn_ppt[3:187]), mean(Wen_05fall_mn_ppt[3:187]), mean(Wen_06fall_mn_ppt[3:187]), mean(Wen_07fall_mn_ppt[3:187]), mean(Wen_08fall_mn_ppt[3:187]), mean(Wen_09fall_mn_ppt[3:187]), mean(Wen_10fall_mn_ppt[3:187]), mean(Wen_11fall_mn_ppt[3:187]), mean(Wen_12fall_mn_ppt[3:187]), mean(Wen_13fall_mn_ppt[3:187]), mean(Wen_14fall_mn_ppt[3:187]), mean(Wen_15fall_mn_ppt[3:187]))
  
  Ann_mn_tmp_sp <- c(mean(Wen_01sp_mn_tmp[3:183]), mean(Wen_02sp_mn_tmp[3:183]), mean(Wen_03sp_mn_tmp[3:183]), mean(Wen_04sp_mn_tmp[3:183]), mean(Wen_05sp_mn_tmp[3:183]), mean(Wen_06sp_mn_tmp[3:183]), mean(Wen_07sp_mn_tmp[3:183]), mean(Wen_08sp_mn_tmp[3:183]), mean(Wen_09sp_mn_tmp[3:183]), mean(Wen_10sp_mn_tmp[3:183]), mean(Wen_11sp_mn_tmp[3:183]), mean(Wen_12sp_mn_tmp[3:183]), mean(Wen_13sp_mn_tmp[3:183]), mean(Wen_14sp_mn_tmp[3:183]), mean(Wen_15sp_mn_tmp[3:183]))
  Ann_mn_tmp_fall <- c(mean(Wen_01fall_mn_tmp[3:187]), mean(Wen_02fall_mn_tmp[3:187]), mean(Wen_03fall_mn_tmp[3:187]), mean(Wen_04fall_mn_tmp[3:187]), mean(Wen_05fall_mn_tmp[3:187]), mean(Wen_06fall_mn_tmp[3:187]), mean(Wen_07fall_mn_tmp[3:187]), mean(Wen_08fall_mn_tmp[3:187]), mean(Wen_09fall_mn_tmp[3:187]), mean(Wen_10fall_mn_tmp[3:187]), mean(Wen_11fall_mn_tmp[3:187]), mean(Wen_12fall_mn_tmp[3:187]), mean(Wen_13fall_mn_tmp[3:187]), mean(Wen_14fall_mn_tmp[3:187]), mean(Wen_15fall_mn_tmp[3:187]))
  
  plot(Ann_mn_tmp_sp, Ann_mn_ppt_sp, pch=c(21,22,24), col="black", bg=col.rainbow, main="Mean sum temp/precip seasonal", xlab = "Mean temp", ylab = "Mean precip", cex=1.5, xlim=c(100,1350), ylim=c(50,600))
  points(Ann_mn_tmp_fall, Ann_mn_ppt_fall, pch=c(21,22,24), col=col.rainbow, bg=col.rainbow, cex=1.5)
  
  
  plot(Ann_mn_tmp_fall, Ann_mn_ppt_fall, pch=c(21,22,24), col="black", bg=col.rainbow, main="Mean sum temp/precip fall", xlab = "Fall mean temp", ylab = "Fallmean precip", cex=1.5, xlim=c(1100,1400), ylim=c(50,250))
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=c(21,22,23), col="black", bg=col.rainbow, bty="n", legend = c("0001", "0102", "0203", "0304", "0405", "0506", "0607", "0708","0809","0910", "1011", "1112", "1213", "1314", "1415"), cex=.7, xpd=NA)
  
  legend("topright", pch=c(15,16,17), col=col.rainbow, bty="n", legend = c("01", "02", "03", "04", "05", "06", "07", "08","09","10", "11", "12", "13", "14", "15"), cex=.7, pt.cex=.95)

#################################################################
#Logger prediction modeling comparison part
  
#################################################################
  
  basin <- "Wen"
  midBasin <- "Wenatchee"
  longBasin <- "Wenatchee"
  dataPath <- "D:/OneDrive/work/research/FishnFire/output/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/Wenatchee/"
  outPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/Climate_analysis/"
  
  coeffs_out <- data.frame(Int=numeric(1), bLST=numeric(1), bLST2=numeric(1), bJul=numeric(1), bElev=numeric(1))
  metrics_out <- data.frame(r2=numeric(1), RMSE=numeric(1), p2=numeric(1), RMSEP=numeric(1), N_Sites=numeric(1), N=numeric(1))
  
  yrPath <- "11"
  yearPath <- "2011"
  #setwd(paste0(dataPath, yearPath, "/Mean"))
  setwd(paste0(mainPath, "/", yearPath))
  
  #Data.in <- read.csv("Model.data.csv", stringsAsFactors = FALSE, header=FALSE)
  Data.in <- read.csv(paste0(basin, "_", yearPath, "_8Day_model_data.csv"), stringsAsFactors = FALSE)
  
  Data.in <- read.csv(paste0(yearPath, "_model_data_Mn.csv"), stringsAsFactors = FALSE)#2011
  
  #colnames(Data.in) <- c("row", "JulDay", "LST", "SiteName", "Date", "Mn", "Mx", "Min", "Year", "GridID", "Elev", "8DMn")
  NoNA.xyz <- Data.in
  
  ind <- apply(Data.in, 1, function(x) !any(is.na(x)))
  NoNA.xyz.all <- Data.in[ind,]
  
  #NoNA.xyz <- NoNA.xyz.all[,c(12, 3, 2, 11, 4)]
  NoNA.xyz <- NoNA.xyz.all[,c(13, 2, 1, 3, 5)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  
  #plot(NoNA.xyz$z, NoNA.xyz$y)
  
  
  y <- NoNA.xyz$y
  x <- NoNA.xyz$x
  z <- NoNA.xyz$z
  e <- NoNA.xyz$e
  #plot(x, y)
  
  
  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  pred.y <- predict(mod)
  #plot(pred.y, y, main = "8-day Mean Full Year")
  #abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  coeffs <- as.matrix(coefficients(mod))
  
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
  
  library(pls)
  mod2 <- plsr(y ~ x + I(x^2) + z + e, validation = "LOO")
  p2 <- R2(mod2)
  detach("package:pls", unload=TRUE)
  
  
  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "Full year"
  pred.out[,5] <- yearPath
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")

  ######
  #pred.out.all <- pred.out
  ######
  pred.out.all <- rbind(pred.out.all, pred.out)
  plot(pred.out.all[,1], pred.out.all[,2])
  abline(0,1)
  
  coeffs_out[1,1] <- coeffs[1,1]
  coeffs_out[1,2] <- coeffs[2,1]
  coeffs_out[1,3] <- coeffs[3,1]
  coeffs_out[1,4] <- coeffs[4,1]
  coeffs_out[1,5] <- coeffs[5,1]
  coeffs_out$Seas <- "Full year"
  coeffs_out$Year <- yearPath
  #####
  #coeffs.out.all <- coeffs_out
  #####
  #coeffs.out.all <- rbind(coeffs.out.all, coeffs_out)
  coeffs.out.all[11,] <- coeffs_out
  
  metrics_out[1,1] <- sum_mod$adj.r.squared
  metrics_out[1,2] <- sum_mod$sigma
  metrics_out[1,3] <- p2$val[5]
  metrics_out[1,4] <- RMSEP
  metrics_out[1,5] <- length(unique(NoNA.xyz$SiteName))
  metrics_out[1,6] <- length(y)
  metrics_out$Seas <- "Full year"
  metrics_out$Year <- yearPath
  #####
  #metrics.out.all <- metrics_out
  #####
  metrics.out.all <- rbind(metrics.out.all, metrics_out)
  
  metrics.out.all[13,] <- metrics_out
  
  
  setwd(outPath)
  write.csv(x=pred.out.all, file=paste0(basin, "_01_15_annual_model_preds.csv"), row.names = FALSE)
  write.csv(x=coeffs.out.all, file=paste0(basin, "_01_15_annual_model_coeffs.csv"), row.names = FALSE)
  write.csv(x=metrics.out.all, file=paste0(basin, "_01_15_annual_model_metrics.csv"), row.names = FALSE)
  
###########################
# water year
############################
  coeffs_out <- data.frame(Int=numeric(1), bLST=numeric(1), bLST2=numeric(1), bJul=numeric(1), bElev=numeric(1))
  metrics_out <- data.frame(r2=numeric(1), RMSE=numeric(1), p2=numeric(1), RMSEP=numeric(1), N_Sites=numeric(1), N=numeric(1))
  
  k <- 1
  year1 <- paste0("200", k)
  yr1 <- paste0("0", k)
  year2 <- paste0("200", k+1)
  yr2 <- paste0("0", k+1)
  
  setwd(paste0(dataPath, year1, "/Mean"))
  Data.in1 <- read.csv("Model.data.csv", stringsAsFactors = FALSE)
  colnames(Data.in1) <- c("row", "JulDay", "LST", "SiteName", "Date", "Mn", "Mx", "Min", "Year", "GridID", "Elev", "8DMn")
  
  setwd(paste0(dataPath, year2, "/Mean"))
  Data.in2 <- read.csv("Model.data.csv", stringsAsFactors = FALSE)
  colnames(Data.in2) <- c("row", "JulDay", "LST", "SiteName", "Date", "Mn", "Mx", "Min", "Year", "GridID", "Elev", "8DMn")
  
  ind <- apply(Data.in1, 1, function(x) !any(is.na(x)))
  NoNA.xyz.all <- Data.in1[ind,]
  NoNA.xyz1 <- NoNA.xyz.all[,c(12, 3, 2, 11, 4)]
  colnames(NoNA.xyz1) <- c("y", "x", "z", "e", "SiteName")
  
  ind <- apply(Data.in2, 1, function(x) !any(is.na(x)))
  NoNA.xyz.all <- Data.in2[ind,]
  NoNA.xyz2 <- NoNA.xyz.all[,c(12, 3, 2, 11, 4)]
  colnames(NoNA.xyz2) <- c("y", "x", "z", "e", "SiteName")
  NoNA.xyz <- rbind(NoNA.xyz1, NoNA.xyz2)
  
  
  