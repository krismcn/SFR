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


  basin <- "JD"
  midBasin <- "JohnDay"
  longBasin <- "JohnDay"
  dataType <- "Temp"
  yrPath <- "00"
  monthPath <- "01"
  yearPath <- "2000"
  dataPath <- "D:/OneDrive/work/GIS/PRISM/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/Climate_analysis/"
    
           
#################################################################
# Data investigation 

# ################################################################
  
  
  setwd(paste0(dataPath, "Tables"))
    
  
  TempMn.out <- NULL
  TempMd.out <- NULL
  Yr.out <- NULL
  PrecipMn.out <- NULL
  PrecipMd.out <- NULL
  
  yrPath <- "09"
  
  
    Temp <- data.frame(mup=NULL)
    for (i in 0:11)
      {
        Temp.in <- read.dbf(paste0(basin, "_Temp_", yrPath, "_", i, ".dbf"))
      Temp <- rbind(Temp, Temp.in)
      }
    
    Temp$Year <- yrPath
    Temp09 <- Temp
    TempMn <- mean(unlist(Temp$Mean)) - 273.15
    TempMd <- median(unlist(Temp$Mean)) - 273.15
    TempMn.out <- append(TempMn.out, TempMn)
    TempMd.out <- append(TempMd.out, TempMd)
    Yr.out <- append(Yr.out, yrPath)
    
    Precip <- data.frame(mup=NULL)
    for (i in 0:11)
    {
      Precip.in <- read.dbf(paste0(basin, "_Precip_", yrPath, "_", i, ".dbf"))
      Precip <- rbind(Precip, Precip.in)
    }
    
    Precip$Year <- yrPath
    Precip09 <- Precip 
    PrecipMn <- mean(unlist(Precip$Mean))
    PrecipMd <- median(unlist(Precip$Mean))
    PrecipMn.out <- append(PrecipMn.out, PrecipMn)
    PrecipMd.out <- append(PrecipMd.out, PrecipMd)
    

##############
# some plots to look at 
# ##############
    
  plot(1:12, Temp00$Mean, main = "Mean Monthly Surface Air Temperature, Wenatchee, 2000-2009", xlab="Month", ylab="Degrees C", ylim=c(265, 295))
  points(1:12, Temp01$Mean, pch=16, col="blue")
  points(1:12, Temp02$Mean, pch=16, col="red")
  points(1:12, Temp03$Mean, pch=16, col="lightblue")
  points(1:12, Temp04$Mean, pch=16, col="orange")
  points(1:12, Temp05$Mean, pch=16, col="pink")
  points(1:12, Temp06$Mean, pch=16, col="purple")
  points(1:12, Temp07$Mean, pch=16, col="palevioletred")
  points(1:12, Temp08$Mean, pch=16, col="grey50")
  points(1:12, Temp09$Mean, pch=16, col="green")
  lines(1:12, Temp00$Mean, col="black")
  lines(1:12, Temp01$Mean, col="blue")
  lines(1:12, Temp02$Mean, col="red")
  lines(1:12, Temp03$Mean, col="lightblue")
  lines(1:12, Temp04$Mean, col="orange")
  lines(1:12, Temp05$Mean, col="pink")
  lines(1:12, Temp06$Mean, col="purple")
  lines(1:12, Temp07$Mean, col="palevioletred")
  lines(1:12, Temp08$Mean, col="grey50")
  lines(1:12, Temp09$Mean, col="green")
  legend("topright", pch=16, title="(K)", legend = c("00","01","02","03","04","05","06","07","08","09"), col=c("black", "blue", "red", "lightblue", "orange", "pink", "purple", "palevioletred", "grey50", "green"), cex=.6)
  
  plot(1:12, Precip00$Mean, main = "Mean Monthly Total Precip, Wenatchee, 2000-2009", xlab="Month", ylab="kgm^2", ylim=c(0, 1.9))
  points(1:12, Precip01$Mean, pch=16, col="blue")
  points(1:12, Precip02$Mean, pch=16, col="red")
  points(1:12, Precip03$Mean, pch=16, col="lightblue")
  points(1:12, Precip04$Mean, pch=16, col="orange")
  points(1:12, Precip05$Mean, pch=16, col="pink")
  points(1:12, Precip06$Mean, pch=16, col="purple")
  points(1:12, Precip07$Mean, pch=21, col="palevioletred")
  points(1:12, Precip08$Mean, pch=16, col="grey50")
  points(1:12, Precip09$Mean, pch=16, col="green")
  lines(1:12, Precip00$Mean, col="black")
  lines(1:12, Precip01$Mean, col="blue")
  lines(1:12, Precip02$Mean, col="red")
  lines(1:12, Precip03$Mean, col="lightblue")
  lines(1:12, Precip04$Mean, col="orange")
  lines(1:12, Precip05$Mean, col="pink")
  lines(1:12, Precip06$Mean, col="purple")
  lines(1:12, Precip07$Mean, col="palevioletred")
  lines(1:12, Precip08$Mean, col="grey50")
  lines(1:12, Precip09$Mean, col="green")
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
        y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, title="kg/m^2", bty="n", legend = c("00","01","02","03","04","05","06","07","08","09"), col=c("black", "blue", "red", "lightblue", "orange", "pink", "purple", "palevioletred", "grey50", "green"), cex=.6, xpd=NA)
  
  mod <- lm(Precip$Mean~Temp$MEAN)
  plot(Temp00$MEAN, Precip$Mean, pch=16, main = "Precipitation by Temperature, Wenatchee, 2000", xlab="Temp", ylab="Precip")
  abline(mod)
  
  plot(Temp00$Mean, Precip00$Mean, main = "Mean Monthly Total Precip by Temp, Wenatchee, 2000-2009", xlab="Month", ylab="kgm^2", ylim=c(0, 1.9), xlim=c(265, 295))
  points(Temp01$Mean, Precip01$Mean, pch=16, col="blue")
  points(Temp02$Mean, Precip02$Mean, pch=16, col="red")
  points(Temp03$Mean, Precip03$Mean, pch=16, col="lightblue")
  points(Temp04$Mean, Precip04$Mean, pch=16, col="orange")
  points(Temp05$Mean, Precip05$Mean, pch=16, col="pink")
  points(Temp06$Mean, Precip06$Mean, pch=16, col="purple")
  points(Temp07$Mean, Precip07$Mean, pch=21, col="palevioletred")
  points(Temp08$Mean, Precip08$Mean, pch=16, col="grey50")
  points(Temp09$Mean, Precip09$Mean, pch=16, col="green")
  abline(lm(Precip00$Mean ~ Temp00$Mean), col="black")
  abline(lm(Precip01$Mean ~ Temp01$Mean), col="blue")
  abline(lm(Precip02$Mean ~ Temp02$Mean), col="red")
  abline(lm(Precip03$Mean ~ Temp03$Mean), col="lightblue")
  abline(lm(Precip04$Mean ~ Temp04$Mean), col="orange")
  abline(lm(Precip05$Mean ~ Temp05$Mean), col="pink")
  abline(lm(Precip06$Mean ~ Temp06$Mean), col="purple")
  abline(lm(Precip07$Mean ~ Temp07$Mean), col="palevioletred")
  abline(lm(Precip08$Mean ~ Temp08$Mean), col="grey50")
  abline(lm(Precip09$Mean ~ Temp09$Mean), col="green")
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00","01","02","03","04","05","06","07","08","09"), col=c("black", "blue", "red", "lightblue", "orange", "pink", "purple", "palevioletred", "grey50", "green"), cex=.6, xpd=NA)
  
  AllTemp <- rbind(Temp00, Temp01, Temp02, Temp03, Temp04, Temp05, Temp06, Temp07, Temp08, Temp09)
  plot(AllTemp$Mean)
  lines(1:120, AllTemp$Mean)
  allPrecip <- rbind(Precip00, Precip01, Precip02, Precip03, Precip04, Precip05, Precip06, Precip07, Precip08, Precip09)
  plot(allPrecip$Mean)
  lines(1:120, allPrecip$Mean)
  
  sumPrecip <- NULL
  sumP <- colSums(Precip00$Mean)
  sumPrecip <- append(sumPrecip, sum(Precip09$Mean))
  plot(TempMn.out, sumPrecip, main="Sum precip by Mean temp")
  
# ###############################
# PRISM raster processing
# ##################################

  #test <- apply(data.out[3:367], 1, cumsum)
  
  basinPts <- read.dbf(paste0(mainPath, basin, "_PRISM_pts.dbf"))
  pts <- data.frame(basinPts[,1], basinPts[,2])
##################
  # precip
# #################
  
  
  for (j in 1998:2015)
    {  
        year <- j
        yr <- j
        setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
        varName <- paste0("Cppt", yr)
        ##### create a list of only the grid files in a directory
        
        allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
        xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
        fileList <- allFiles[!allFiles %in% xmlFiles]
        
        ###### read in one grid to get the structure
        
        #r <- readGDAL(fileList[1])
        rRaster <- raster(fileList[1])
        
        ##### clip the raster to a reasonable PNW extent
        
        # ext <- extent(-125, -107, 40, 50)
        # rExt <- crop(rRaster, ext)
        # 
        ##### convert the raster to points and build the data structure for the loop
        
        # rPoints <- rasterToPoints(rExt)
        # 
       
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
        
        data.out$PtID <- basinPts$PtID
        write.dbf(data.out, file= paste0("SumPpt", year, ".dbf"))
        
        
        assign(paste0("Ppt_", basin, "_", year), colMeans(data.out))
        
    }

  list <- mget(paste0("Ppt_", basin, "_", 1981:2010))
  all <- do.call(rbind, list)
  
    
  ############# that one time I needed to reread everything ##########
  
  for (j in 2000:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
      test.in <- read.dbf(paste0("SumPpt", year, ".dbf"))
      data.out_JD <- test.in[test.in$PtID %in% basinPts$PtID,]
      assign(paste0("Ppt_JD_", year), colMeans(data.out_JD))
    }
  #################### This part just creates a list of all the vector names and row binds them ##############
  
  list <- mget(paste0("Ppt_JD_", 2000:2015))
  all <- do.call(rbind, list)
  
  y01_15_ppt <- all[2:16, 367]
  write.csv(y01_15_ppt, file=paste0(basin, "_01_15_mn_sum_precip.csv"))
######################
# if going back later to process
# #####################
  
  for (j in 1981:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
      data.in <- read.dbf(paste0("SumTmp", year, ".dbf"))
      data.out <- data.in[data.in$PtID %in% basinPts$PtID,]
      write.dbf(data.out, file= paste0("SumTmpJD", year, ".dbf"))
    }
    
  
  
  y30_Ppt_mn <- colMeans(all[,3:367])
  
  
# #####
# # Subsetting and summarizing by basin ppt
# #####
#   
  y30_precip <- c(mean(Ppt_JD_1981), mean(Ppt_JD_1982), mean(Ppt_JD_1983), mean(Ppt_JD_1984), mean(Ppt_JD_1985), mean(Ppt_JD_1986), mean(Ppt_JD_1987), mean(Ppt_JD_1988), mean(Ppt_JD_1989), mean(Ppt_JD_1990), mean(Ppt_JD_1991), mean(Ppt_JD_1992), mean(Ppt_JD_1993), mean(Ppt_JD_1994), mean(Ppt_JD_1995), mean(Ppt_JD_1996), mean(Ppt_JD_1997), mean(Ppt_JD_1998), mean(Ppt_JD_1999), mean(Ppt_JD_2000), mean(Ppt_JD_2001), mean(Ppt_JD_2002), mean(Ppt_JD_2003), mean(Ppt_JD_2004), mean(Ppt_JD_2005), mean(Ppt_JD_2006), mean(Ppt_JD_2007), mean(Ppt_JD_2008), mean(Ppt_JD_2009), mean(Ppt_JD_2010), mean(Ppt_JD_2011), mean(Ppt_JD_2012), mean(Ppt_JD_2013), mean(Ppt_JD_2014), mean(Ppt_JD_2015))
  
  
  y01_15_precip <- c(mean(Ppt_JD_2001), mean(Ppt_JD_2002), mean(Ppt_JD_2003), mean(Ppt_JD_2004), mean(Ppt_JD_2005), mean(Ppt_JD_2006), mean(Ppt_JD_2007), mean(Ppt_JD_2008), mean(Ppt_JD_2009), mean(Ppt_JD_2010), mean(Ppt_JD_2011), mean(Ppt_JD_2012), mean(Ppt_JD_2013), mean(Ppt_JD_2014), mean(Ppt_JD_2015))
  
  list <- mget(paste0("Ppt_", basin, "_", 1981:2010))
  all <- do.call(rbind, list)
  y30_Ppt_mn <- colMeans(all[,3:367])
  
  
  plot(Ppt_JD_2000[3:367], pch=16, col="cornflowerblue", cex=.7, main = "Mean cumulative precip, John Day basin, 2000-2015", xlab="Julian Day", ylab="Precip", ylim=c(0,600))
  points(Ppt_JD_2001[3:367], pch=16, col="cyan4", cex=.7)
  points(Ppt_JD_2002[3:367], pch=16, col="khaki", cex=.7)
  points(Ppt_JD_2003[3:367], pch=16, col="deeppink4", cex=.7)
  points(Ppt_JD_2004[3:367], pch=16, col="deeppink", cex=.7)
  points(Ppt_JD_2009[3:367], pch=16, col="darkgrey", cex=.7)
  points(Ppt_JD_2008[3:367], pch=16, col="lightblue", cex=.7)
  points(Ppt_JD_2007[3:367], pch=16, col="red", cex=.7)
  points(Ppt_JD_2011[3:367], pch=16, col="palevioletred", cex=.7)
  points(Ppt_JD_2006[3:367], pch=16, col="darkgoldenrod1", cex=.7)
  points(Ppt_JD_2005[3:367], pch=16, col="blueviolet", cex=.7)
  points(Ppt_JD_2010[3:367], pch=16, col="lightgrey", cex=.7)
  points(Ppt_JD_2012[3:367], pch=16, col="darkolivegreen4", cex=.7)
  points(Ppt_JD_2013[3:367], pch=16, col="darkturquoise", cex=.7)
  points(Ppt_JD_2014[3:367], pch=16, col="chartreuse2", cex=.7)
  points(Ppt_JD_2015[3:367], pch=15, col="blue2", cex=.7)

  points(y30_Ppt_mn[1:365], pch=16, col="black", cex=.8)
  
   legend(x=grconvertX(c(1.0,1.4), from='npc'),
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)

  
##################
# Temperature
# #################
  
  
  for (j in 2010:2011)
    {  
      year <- j
      yr <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
      varName <- paste0("Ctmp", yr)
      ##### create a list of only the grid files in a directory
      
      allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
      xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
      fileList <- allFiles[!allFiles %in% xmlFiles]
      
      ###### read in one grid to get the structure
      
      rRaster <- raster(fileList[1])
      
    
      
      data.out <- data.frame(pts, extract(rRaster, pts))
      colnames(data.out) <- c("x", "y", paste0(varName,"001"))
      data.out[,3][data.out[,3] < 0]  <- 0
      
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
      
      data.out$PtID <- basinPts$PtID
      write.dbf(data.out, file= paste0("SumTmp", year, ".dbf"))
      
      
      assign(paste0("Tmp_", basin, "_", year), colMeans(data.out))
      
    }

  
  y30_temp <- c(mean(Tmp_JD_1981), mean(Tmp_JD_1982), mean(Tmp_JD_1983), mean(Tmp_JD_1984), mean(Tmp_JD_1985), mean(Tmp_JD_1986), mean(Tmp_JD_1987), mean(Tmp_JD_1988), mean(Tmp_JD_1989), mean(Tmp_JD_1990), mean(Tmp_JD_1991), mean(Tmp_JD_1992), mean(Tmp_JD_1993), mean(Tmp_JD_1994), mean(Tmp_JD_1995), mean(Tmp_JD_1996), mean(Tmp_JD_1997), mean(Tmp_JD_1998), mean(Tmp_JD_1999), mean(Tmp_JD_2000), mean(Tmp_JD_2001), mean(Tmp_JD_2002), mean(Tmp_JD_2003), mean(Tmp_JD_2004), mean(Tmp_JD_2005), mean(Tmp_JD_2006), mean(Tmp_JD_2007), mean(Tmp_JD_2008), mean(Tmp_JD_2009), mean(Tmp_JD_2010), mean(Tmp_JD_2011), mean(Tmp_JD_2012), mean(Tmp_JD_2013), mean(Tmp_JD_2014), mean(Tmp_JD_2015))
  
  y01_15_temp <- c(mean(Tmp_JD_2001), mean(Tmp_JD_2002), mean(Tmp_JD_2003), mean(Tmp_JD_2004), mean(Tmp_JD_2005), mean(Tmp_JD_2006), mean(Tmp_JD_2007), mean(Tmp_JD_2008), mean(Tmp_JD_2009), mean(Tmp_JD_2010), mean(Tmp_JD_2011), mean(Tmp_JD_2012), mean(Tmp_JD_2013), mean(Tmp_JD_2014), mean(Tmp_JD_2015))
  
  list <- mget(paste0("Tmp_", basin, "_", 1981:2010))
  all <- do.call(rbind, list)
  y30_Tmp <- colMeans(all[,3:367])
  
  plot(Tmp_JD_2000[3:367], pch=16, col="cornflowerblue", cex=.7, main = "Mean cumulative temp, John Day basin, 2000-2015", xlab="Julian Day", ylab="Precip", ylim=c(0,4000))
  points(Tmp_JD_2001[3:367], pch=16, col="cyan4", cex=.7)
  points(Tmp_JD_2002[3:367], pch=16, col="khaki", cex=.7)
  points(Tmp_JD_2003[3:367], pch=16, col="deeppink4", cex=.7)
  points(Tmp_JD_2004[3:367], pch=16, col="deeppink", cex=.7)
  points(Tmp_JD_2009[3:367], pch=16, col="darkgrey", cex=.7)
  points(Tmp_JD_2008[3:367], pch=16, col="lightblue", cex=.7)
  points(Tmp_JD_2007[3:367], pch=16, col="red", cex=.7)
  points(Tmp_JD_2011[3:367], pch=16, col="palevioletred", cex=.7)
  points(Tmp_JD_2006[3:367], pch=16, col="darkgoldenrod1", cex=.7)
  points(Tmp_JD_2005[3:367], pch=16, col="blueviolet", cex=.7)
  points(Tmp_JD_2010[3:367], pch=16, col="lightgrey", cex=.7)
  points(Tmp_JD_2012[3:367], pch=16, col="darkolivegreen4", cex=.7)
  points(Tmp_JD_2013[3:367], pch=16, col="darkturquoise", cex=.7)
  points(Tmp_JD_2014[3:367], pch=16, col="chartreuse2", cex=.7)
  points(Tmp_JD_2015[3:367], pch=15, col="blue2", cex=.7)
  
  points(y30_Tmp[3:367], pch=16, col="black", cex=.8)
  
  legend(x=grconvertX(c(1.0,1.4), from='npc'),
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("00", "01", "02", "03", "04", "05", "06", "07","08","09", "10", "11", "12", "13", "14", "15"), col=c("cornflowerblue", "cyan4", "khaki", "deeppink4", "deeppink", "blueviolet", "darkgoldenrod1", "red", "lightblue", "black", "lightgrey", "palevioletred", "darkolivegreen4", "darkturquoise", "chartreuse2", "blue2"), cex=.65, xpd=NA)
  
 
############# that one time I needed to reread everything ##########
  
  for (j in 1982:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
      data.in <- read.dbf(paste0("SumTmp", year, ".dbf"))
      data.in[data.in < 0] <- 0
      assign(paste0("Tmp_", basin, "_", year), colMeans(data.in))
  }
  
############# that one time I needed to rename everything ##########
  
  for (j in 1981:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
      data.in <- read.dbf(paste0("SumPpt", year, ".dbf"))
      write.dbf(data.out, file= paste0("SumPptJD", year, ".dbf"))
    }
#################### This part just creates a list of all the vector names and row binds them ##############
  
  lister <- mget(paste0("Ppt_JD_", 1981:2015))
  all <- do.call(rbind, list)
  y30_precip <- colMeans(all)
  
  y30_Ppt_mn<- mean(y30_precip[3:367])
  col.rainbow <- rainbow(15)
  
  plot(y01_15_temp, y01_15_precip, pch=c(21,22,24), col="black", bg=col.rainbow, main="Mean sum temp/precip (calendar year)", xlab = "Mean temp", ylab = "Mean precip", cex=1.5,)
  points(y30_Tmp_mn, y30_Ppt_mn, pch=8, col="black", cex=1.5)
  

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
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/Temp/", year))
      data.in <- read.dbf(paste0("SumTmp", basin, year, ".dbf"))
      assign(paste0("Tmp_", basin, "_", year), colMeans(data.in))
    }
  
  list <- mget(paste0("Tmp_JD_", 1981:2010))
  all <- do.call(rbind, list)
  y30_Tmp <- colMeans(all[,3:367])
  
  y30_mn_temp <- mean(y30_Tmp)
  y30_mn_precip <- mean(y30_Ppt)
  
  i <- 1
  HUCpts <- basinPts[basinPts$HUC == i, "PtID"]
  
  for (j in 1981:2015)
    {
      year <- j
      setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
      data.in <- read.dbf(paste0("SumPpt", basin, year, ".dbf"))
      data.out <- data.in[data.in$PtID %in% HUCpts,]
      assign(paste0("Ppt_", basin, "_HUC_", i, "_", year), colMeans(data.out))
    }
  
    
 
  lister <- mget(paste0("Tmp_", basin, "_HUC_", i, "_", 1981:2010))
  all <- do.call(rbind, lister)
  assign(paste0("y30_Tmp_HUC", i), colMeans(all))
  
  y30_Tmp_HUC3_mn <- mean(y30_Tmp_HUC3[3:367])
  y30_Ppt_HUC3_mn <- mean(y30_Ppt_HUC3[3:367])
  
  Ann_mn_ppt_HUC1 <- c(mean(Ppt_JD_HUC_1_2000[3:367]), mean(Ppt_JD_HUC_1_2001[3:367]), mean(Ppt_JD_HUC_1_2002[3:367]), mean(Ppt_JD_HUC_1_2003[3:367]), mean(Ppt_JD_HUC_1_2004[3:367]), mean(Ppt_JD_HUC_1_2005[3:367]), mean(Ppt_JD_HUC_1_2006[3:367]), mean(Ppt_JD_HUC_1_2007[3:367]), mean(Ppt_JD_HUC_1_2008[3:367]), mean(Ppt_JD_HUC_1_2009[3:367]), mean(Ppt_JD_HUC_1_2010[3:367]), mean(Ppt_JD_HUC_1_2011[3:367]), mean(Ppt_JD_HUC_1_2012[3:367]), mean(Ppt_JD_HUC_1_2013[3:367]), mean(Ppt_JD_HUC_1_2014[3:367]), mean(Ppt_JD_HUC_1_2015[3:367]))
  
  Ann_mn_tmp_HUC1<- c(mean(Tmp_JD_HUC_1_2000[3:367]), mean(Tmp_JD_HUC_1_2001[3:367]), mean(Tmp_JD_HUC_1_2002[3:367]), mean(Tmp_JD_HUC_1_2003[3:367]), mean(Tmp_JD_HUC_1_2004[3:367]), mean(Tmp_JD_HUC_1_2005[3:367]), mean(Tmp_JD_HUC_1_2006[3:367]), mean(Tmp_JD_HUC_1_2007[3:367]), mean(Tmp_JD_HUC_1_2008[3:367]), mean(Tmp_JD_HUC_1_2009[3:367]), mean(Tmp_JD_HUC_1_2010[3:367]), mean(Tmp_JD_HUC_1_2011[3:367]), mean(Tmp_JD_HUC_1_2012[3:367]), mean(Tmp_JD_HUC_1_2013[3:367]), mean(Tmp_JD_HUC_1_2014[3:367]), mean(Tmp_JD_HUC_1_2015[3:367]))
 
  
  Ann_mn_HUC_1 <- data.frame(Ppt=Ann_mn_ppt_HUC1, Tmp=Ann_mn_tmp_HUC1, year=sprintf('%02d', 0:15), HUC=1)
  
  plot(Ann_mn_HUC_1$Tmp, Ann_mn_HUC_1$Ppt, pch=21, col="black", bg="blue", main="Mean sum temp/precip JD (by HUC) 2000-2015", xlab = "Mean temp", ylab = "Mean precip", cex=1.5, xlim=c(1100,2100), ylim=c(100,350),)
  points(Ann_mn_HUC_2$Tmp, Ann_mn_HUC_2$Ppt, pch=22, col="black", bg="red", cex=1.5)
  points(Ann_mn_HUC_3$Tmp, Ann_mn_HUC_3$Ppt, pch=24, col="black", bg="green", cex=1.5)
  points(Ann_mn_HUC_4$Tmp, Ann_mn_HUC_4$Ppt, pch=23, col="black", bg="yellow", cex=1.5)
  points(y30_Tmp_HUC1_mn, y30_Ppt_HUC1_mn, pch="#", col="blue", cex=1.5)
  points(y30_Tmp_HUC2_mn, y30_Ppt_HUC2_mn, pch="#", col="red", cex=1.5)
  points(y30_Tmp_HUC3_mn, y30_Ppt_HUC3_mn, pch="#", col="green", cex=1.5)
  points(y30_Tmp_HUC4_mn, y30_Ppt_HUC4_mn, pch="#", col="yellow", cex=1.5)
  points(y30_mn_temp, y30_mn_precip, pch=8, col="black", cex=2.5)
  
  legend("topright", pch=c(16,15,17,18,8), col=c("blue", "red", "green", "yellow", "black"), bg=c("blue", "red", "green", "yellow"),bty="n", legend = c("SFJD", "NFJD", "MFJD", "LFJD", "30yJD"), cex=.7, pt.cex=.95)

  #Ann_mn_HUC_4$year[1:15] <- ""  
  
  textxy(Ann_mn_HUC_1$Tmp, Ann_mn_HUC_1$Ppt, Ann_mn_HUC_1$year, cex=.8, offset = 1.0, pos=1)
  
  
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
  
  
  