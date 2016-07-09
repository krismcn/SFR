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

#################################################################
  
  
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
###############
    
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
  
################################
# PRISM raster processing
###################################

  
  library(rgdal)
  library(raster)
  library(rNOMADS)
  library(gdalUtils)

  year <- "2011"
  setwd(paste0("D:/OneDrive/work/GIS/PRISM/", year))
  varName <- "Cppt11"
##### create a list of only the grid files in a directory
  
  allFiles <- list.files(pattern="*bil.bil", full.names=TRUE)
  xmlFiles <- list.files(pattern="*bil.bil.aux.xml", full.names=TRUE)
  fileList <- allFiles[!allFiles %in% xmlFiles]
  
###### read in one grid to get the structure
  
  r <- readGDAL(fileList[1])
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
  write.dbf(data.out, "SumPpt2011.dbf")
  data.out11 <- data.out
  
#####
# Subsetting and summarizing by basin
#####
  
  basinPts <- read.dbf("D:/OneDrive/work/GIS/PRISM/Wen_prism_pts.dbf")
  data.out11_Wen <- data.out11[data.out11$PtID %in% basinPts$PtID,]
  Wen_11_mn <- colMeans(data.out11_Wen)
  plot(Wen_10_mn[3:367], pch=16, col="lightgrey", main = "Mean cumulative precip (mm), Wenatchee basin, 2007-2010", xlab="Julian Day", ylab="Precip")
  points(Wen_09_mn[3:367], pch=16, col="black")
  points(Wen_08_mn[3:367], pch=16, col="lightblue")
  points(Wen_07_mn[3:367], pch=16, col="red")
  points(Wen_11_mn[3:367], pch=16, col="palevioletred")
  legend(x=grconvertX(c(1.0,1.4), from='npc'), 
         y=grconvertY(c(0.6, 0.8), from='npc'), pch=16, bty="n", legend = c("07","08","09", "10", "11"), col=c("red", "lightblue", "black", "lightgrey", "palevioletred"), cex=.6, xpd=NA)
  
###########################################
    
  
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5
  
  pred.new <- predict(mod, data = data.sp)
  pred.new[pred.new < -0.5] = -0.5
  
  data.sp$pred <- unlist(pred.new)
  plot(pred.y, y, main = "8-day Mn Spring Leg")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  gvmodel <- gvlma(mod)
  summary(gvmodel)
  plot(mod, which= 1:5)
  outlierTest(mod)
  qqPlot(mod, main="QQ Plot Spring Leg")
  spreadLevelPlot(mod)
  plot(pred.y, mod$residuals, main="Model diagnostics Spring Leg", xlab="Predicted", ylab="Residuals")
  
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
  pred.out[,4] <- "Spring"
  pred.out[,5] <- yearPath
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_Mean_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  
  
  coeffs_out[1,1] <- coeffs[1,1]
  coeffs_out[1,2] <- coeffs[2,1]
  coeffs_out[1,3] <- coeffs[3,1]
  coeffs_out[1,4] <- coeffs[4,1]
  coeffs_out[1,5] <- coeffs[5,1]
  
  metrics_out[1,1] <- sum_mod$adj.r.squared
  metrics_out[1,2] <- sum_mod$sigma
  metrics_out[1,3] <- p2$val[5]
  metrics_out[1,4] <- RMSEP
  metrics_out[1,5] <- length(unique(data.sp$SiteName))
  metrics_out[1,6] <- length(y)



  y <- data.fall$y
  x <- data.fall$x
  z <- data.fall$z
  e <- data.fall$e
  plot(z, y)  
  plot(x, y, main="fall")

  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  sum_mod
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  pressstat_sum$stat


###########################################

  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5
  
  pred.new <- predict(mod, data=data.fall)
  pred.new[pred.new < -0.5] = -0.5
  
  data.fall$pred <- unlist(pred.new)
  
  plot(pred.y, y, main = "8-day Mean stream temperature (Fall Leg)", xlab="Predicted ('C)", ylab="Observed ('C)")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  
  plot(mod, which = 1:5)
  gvmodel <- gvlma(mod)
  summary(gvmodel)
  outlierTest(mod)
  qqPlot(mod, main="QQ Plot Fall Leg")
  spreadLevelPlot(mod)
  plot(pred.y, mod$residuals, main="Model diagnostics Fall Leg", xlab="Predicted", ylab="Residuals")
  
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
  pred.out[,4] <- "fall"
  pred.out[,5] <- yearPath
  plot(pred.out[,1], pred.out[,2])
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_Max_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  

  coeffs_out[2,1] <- coeffs[1,1]
  coeffs_out[2,2] <- coeffs[2,1]
  coeffs_out[2,3] <- coeffs[3,1]
  coeffs_out[2,4] <- coeffs[4,1]
  coeffs_out[2,5] <- coeffs[5,1]

  metrics_out[2,1] <- sum_mod$adj.r.squared
  metrics_out[2,2] <- post_mod$sigma
  metrics_out[2,3] <- p2$val[5]
  metrics_out[2,4] <- RMSEP
  metrics_out[2,5] <- length(unique(data.fall$SiteName))
  metrics_out[2,6] <- length(y)

write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mn.csv"), sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_Mn.csv"), sep = ",", col.names=T)

  pred.y <- read.csv(paste0("jk_pred_v_y_Mean_", basin, "_", yearPath, "_sp_fall.csv"), stringsAsFactors = FALSE)
  colnames(pred.y) <- c("Y", "PredY", "JulDay", "Season", "Year")
  
  plot(pred.y$PredY, pred.y$Y, pch=16, col="blue", main=paste0("Mn 8-day stream temp ", basin, " ", yearPath), xlab="Predicted", ylab="Observed")
  abline(0,1)
  abline(lm(pred.y$Y~ pred.y$PredY), col="blue")
  fit <- lm(pred.y$Y~ pred.y$PredY)
  plot(fit)
  summary(fit)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################
setwd(paste0(mainPath, longBasin))  

  elev.in <- read.csv(paste0(basin, "_rca_elev.csv"), stringsAsFactors = F)

setwd(paste0(mainPath, longBasin, "/"))
  
  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))

  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])

  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")
  

setwd(paste0(mainPath, longBasin, "/", yearPath))

  coeffs.in <- read.csv(paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mn.csv"), stringsAsFactors=FALSE)
  LogPred.out <- LST.elev[,c(2:47,48,1)]
  LogPred.out[,1:46] <- 0
  rcas <- unique(LST.elev$RCAID)
  LST.sum <- LST.elev[,c(2:47, 48,1)]

  for (i in 1:length(rcas))  
    {
      x <- unlist(LST.sum[i,])
      maxrow <- as.numeric(which.max(x[1:46])) #either specify or let be dynamic
      midrow <- maxrow - 1
      day <- as.numeric(colnames(LST.sum)[maxrow])
      
      j <- 1
      for (l in 1:midrow)
      {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[47] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
       j <- j + 8}
      k <- day
      for (l in maxrow:46)     
      {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[47] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
       k <- k + 8}
      if (maxrow > 3)
        {
        x[(midrow)] <- NA
        fill <- na.spline(x[1:46],)
        x[(maxrow-3):(maxrow+3)] <- fill[(maxrow-3):(maxrow+3)]
        }
      LogPred.out[i,1:46] <- x [1:46] 
    }


  LogPred.out <- as.data.frame(LogPred.out)

  LogPred.out$Basin_RCA <- paste0(basin, "_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:46]))
  varName <- paste0("Tmn_", yrPath)
  names.out <- sprintf("%s_%03d", varName, namesnum)
  colnames(LogPred.out)[1:46] <- names.out[1:46]
  plot(namesnum, LogPred.out[1,1:46])

  LogPred.out[LogPred.out < -0.5] = -0.5

write.dbf(LogPred.out, file = paste0("predt", yearPath, "_", basin, "_8D_Mn.dbf")) 

##########################
#This parts formats the error by day/site info
###########################

  # NoNA.xyz <- read.csv(paste0(basin, "_", yearPath, "_8Day_model_data.csv"), stringsAsFactors = TRUE)
  # NoNA.xyz <- orderBy(~z, NoNA.xyz)
  # maxrow <- which.max(NoNA.xyz$y)
  # data.sp <- NoNA.xyz[1:(maxrow-1),]
  # data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]
  # 
  # y <- data.sp$y
  # x <- data.sp$x
  # z <- data.sp$z
  # e <- data.sp$e
  # mod <- lm(y ~ x + I(x^2) + z + e)
  # pred.new <- predict(mod, data = data.sp)
  # pred.new[pred.new < -0.5] = -0.5
  # data.sp$pred <- unlist(pred.new)
  # 
  # y <- data.fall$y
  # x <- data.fall$x
  # z <- data.fall$z
  # e <- data.fall$e
  # mod <- lm(y ~ x + I(x^2) + z + e)
  # pred.new <- predict(mod, data=data.fall)
  # pred.new[pred.new < -0.5] = -0.5
  # data.fall$pred <- unlist(pred.new)
  # 
  error.pts <- rbind(data.sp, data.fall)
  SiteID <- unique(error.pts$SiteName)
  SiteID <- as.matrix(SiteID)
  Error.pts.out <- matrix(nrow = length(SiteID), ncol = 47)
  
  errorName <- paste0("JulDay_", namesnum)
  colnames(Error.pts.out)[2:47] <- errorName
  colnames(Error.pts.out)[1] <- "SiteName"
  
  Error.pts.out[1:length(SiteID), 1] <- unlist(SiteID)[1:length(SiteID)]
  Error.pts.out <- as.data.frame(Error.pts.out, stringsAsFactors = FALSE)


  for (i in SiteID) 
    { 
      error.site <- error.pts[error.pts$SiteName  == i,]
      error.site$error <- error.site[,'y']-error.site[,'pred']
      
      
      error <- matrix(ncol=1, nrow=46)
      error[,1] <- namesnum
      error <- data.frame(error)
      colnames(error) <- c("JulDay")
      
      error.site.fill <- merge(error, error.site, by.x = "JulDay", by.y = "z", all.x=TRUE, all.y = FALSE)
      Error.pts.out[Error.pts.out$SiteName==i,2:47] <- as.numeric(unlist(error.site.fill$error))
    }

Error.pts.out[,2:47] <- sapply(Error.pts.out[,2:47], as.numeric)

Error.pts.out[,2:47] <- round(Error.pts.out[,2:47], digits=3)
plot(1:46, Error.pts.out[2,2:47])
write.dbf(Error.pts.out, file = paste0("Error", yearPath, "_", basin, "_8D_Mn.dbf")) 
write.csv(Error.pts.out, file = paste0("Error", yearPath, "_", basin, "_8D_Mn.csv"))

###############################################
#
# Edited 5 may 2016 to add animation output
# Edited 25 May 2016 to add HUC display
################################################

library(rgdal)
library(RColorBrewer)
library(classInt)
library(TeachingDemos)

  modelPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", yearPath, "_temp_CHaMP/")
  ptsPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", "Mean_models")
  shpPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yearPath, "/", "graphic_shapes")

  netname <- paste0(basin, "_", yearPath, "_8D_mn")
  ptsname <- paste0(basin, "_Error_", yearPath, "_8D_mn")
  
setwd(ptsPath)

  error_pts <- readOGR(dsn=".", "Wen_Error_2014_8D_Mn")
  
setwd(shpPath) 
  
  network <- readOGR(dsn=".", layer = netname)
  network@data <- network@data[,-3]
  
  HUCs <- readOGR(dsn=".", "Wen_HUC5")
  
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)

  names.out <- colnames(network@data[4:49])
  namesnum <- as.numeric(gsub("Tmn_14_", "", colnames(network@data[4:49])))
  means <- colMeans(network@data[4:49])
  SDs <- colStdevs(network@data[4:49])
  yplus <- means + SDs
  yminus <- means - SDs
  df <- data.frame(means=means, SDs=SDs, names=namesnum)
  sequ <- c(1:46)
  namer <- sprintf('%03d', sequ)
  fix4 <- classIntervals(means, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix4.colors <- findColours(fix4,pal=seis)
  
  for (i in 4:49)
    {
      
      namey <- gsub("Tmn_14_", "", colnames(network@data)[i])
      
      filename <- paste0(mainPath, longBasin, "/", yearPath, "/graphics3/", namer[i-3], ".png", sep="")
      png(filename=filename, res = 300, width = 1500, height = 1500, units = "px", bg="black")
      
      fix3 <- classIntervals(network@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix3.colors <- findColours(fix3,pal=seis)
      
      cexEr <- ifelse(abs(error_pts@data[,i]) <= 1, 0.5,
                      ifelse(abs(error_pts@data[,i])>1, 0.75,
                             ifelse(abs(error_pts@data[,i])>2, 1.0,
                                    ifelse(abs(error_pts@data[,i])>3, 1.25, NA))))
      
      plot(network, col=fix3.colors, bg="black", fg="white")
      points(error_pts, pch=16, col="gray40", cex=cexEr)
      
      legend("right", fill = attr(fix3.colors, "palette"), legend = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18+"), bty = "n", cex=.5, inset=c(.1,0), text.col="white");
      
      
      title("Wenatchee 8-day mean 2014 ('C)", sub = paste0("Julian Day ", namey), line=-0.9, adj=.80, col.main="white", col.sub="white", outer=FALSE, cex.main=0.5, cex.sub=0.5)
      tmp2 <- subplot(
        plot(namesnum[1:(i-2)], means[1:(i-2)], col=fix4.colors, pch=16, bty="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean (+/-SD)", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"), 
        x=grconvertX(c(0.1,0.45), from='npc'), 
        y=grconvertY(c(0.05, 0.20), from='npc'),
        size=c(1,1.5), vadj=0.5, hadj=0.5, 
        pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
      
      
      dev.off()
    }


setwd(paste0(mainPath, longBasin, "/", yearPath, "/graphics3/"))

system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png example2.mpeg')

###############################################
#
# Edited 5 may 2016 to add animation output
# Edited 25 May 2016 to add HUC display
# Edited 13 June 2016 to add multiple frame display
################################################

library(rgdal)
library(RColorBrewer)
library(classInt)

  yr1 <- 2012
  yr2 <- 2013
  yr3 <- 2014
  
  ptsPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", "Mean_models")
  shpPath3 <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yr3, "/", "graphic_shapes")
  shpPath2 <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yr2, "/", "graphic_shapes")
  
  netname <- paste0(basin, "_base_net")
  netname1 <- paste0(basin, "_", yr1, "_8D_mn")
  netname2 <- paste0(basin, "_", yr2, "_8D_mn") 
  netname3 <- paste0(basin, "_", yr3, "_8D_mn")
  
  ptsname1 <- paste0(basin, "_Error_", yr1, "_8D_Mn")
  ptsname2 <- paste0(basin, "_Error_", yr2, "_8D_Mn")
  ptsname3 <- paste0(basin, "_Error_", yr3, "_8D_Mn")

setwd(ptsPath)

  error_pts1 <- readOGR(dsn=".", layer = ptsname1)
  
  error_pts3 <- readOGR(dsn=".", layer = ptsname3)

setwd(shpPath1) 

  network <- readOGR(dsn=".", layer = netname)
  
  network1 <- readOGR(dsn=".", layer = netname1)
  network1@data <- network1@data[,-4]

setwd(shpPath2)

  network2 <- readOGR(dsn=".", layer = netname2)
  network2@data <- network2@data[,-4]
  error_pts2 <- readOGR(dsn=".", layer = ptsname2)
  
setwd(shpPath3)

  network3 <- readOGR(dsn=".", layer = netname3)
  
  
  #HUCs <- readOGR(dsn=".", "Wen_HUC5")
  
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)
  
  ##################### basin means ####################
  
  means3 <- colMeans(network3@data[4:49])
  SDs3 <- colStdevs(network3@data[4:49])
  yplus3 <- means3 + SDs3
  yminus3 <- means3 - SDs3
  df3 <- data.frame(means=means3, SDs=SDs3, names=namesnum)
  
  fix4 <- classIntervals(means3, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix4.colors <- findColours(fix4,pal=seis)
  
  plot(namesnum[1:(i-2)], means3[1:(i-2)], col=fix4.colors, pch=16, bty="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean (+/-SD)", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black")
  
  means2 <- colMeans(network2@data[4:49])
  SDs2 <- colStdevs(network2@data[4:49])
  yplus2 <- means2 + SDs2
  yminus2 <- means2 - SDs2
  df2 <- data.frame(means=means2, SDs=SDs2, names=namesnum)
  points(namesnum[1:(i-2)], means2[1:(i-2)], col=fix4.colors, pch=15, bty="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean (+/-SD)", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black")
  
  ######################################################
  
  sequ <- c(1:46)
  namer <- sprintf('%03d', sequ)
  
  
  names.out <- colnames(network3@data[4:49])
  namesnum <- as.numeric(gsub("Tmn_14_", "", colnames(network3@data[4:49])))

  op <- par(mfrow = c(1,3),
            mar = c(0,0,1,1) + 0.1)
  usr <- par("usr")
  
  for (i in 25:26)
    {
      
      namey <- gsub("Tmn_14_", "", colnames(network3@data)[i])
      
      filename <- paste0(mainPath, longBasin, "/", yr3, "/graphics5/", namer[i-3], ".png", sep="")
      png(filename=filename, res = 300, width = 1500, height = 1500, units = "px", bg="black")
      
      fix3 <- classIntervals(network3@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix3.colors <- findColours(fix3,pal=seis)
      
      # cexEr1 <- ifelse(abs(error_pts1@data[,i-1]) <= 1, 0.5,
      #                 ifelse(abs(error_pts1@data[,i-1])>1, 0.75,
      #                        ifelse(abs(error_pts1@data[,i-1])>2, 1.0,
      #                               ifelse(abs(error_pts1@data[,i-1])>3, 1.25, NA))))
      
      cexEr2 <- ifelse(abs(error_pts2@data[,i-1]) <= 1, 0.5,
                       ifelse(abs(error_pts2@data[,i-1])>1, 0.75,
                              ifelse(abs(error_pts2@data[,i-1])>2, 1.0,
                                     ifelse(abs(error_pts2@data[,i-1])>3, 1.25, NA))))
      
      cexEr3 <- ifelse(abs(error_pts3@data[,i-1]) <= 1, 0.5,
                       ifelse(abs(error_pts3@data[,i-1])>1, 0.75,
                              ifelse(abs(error_pts3@data[,i-1])>2, 1.0,
                                     ifelse(abs(error_pts3@data[,i-1])>3, 1.25, NA))))
      
      fix2 <- classIntervals(network2@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix2.colors <- findColours(fix2,pal=seis)
      
      fix1 <- classIntervals(network1@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix1.colors <- findColours(fix1,pal=seis)
      
      op <- par(mfrow = c(1,3),
                oma = c(0,0,0,4),
                mar = c(0,0,0,0) + 0.1,
                mgp = c(2,0,0),
                xpd=NA)
      
     
      plot(network, col="gray25", bg="black")
      plot(network1, col=fix1.colors, bg="black", fg="white", add=TRUE)
      #points(error_pts1, pch=16, col="gray40", cex=cexEr1)
      text(x = (usr[1] + usr[2])/2, y=usr[3], yr1, col="white", cex=.8)
      
      plot(network, col="gray25", bg="black")
      plot(network2, col=fix2.colors, bg="black", fg="white", add=TRUE)
      points(error_pts2, pch=16, col="gray40", cex=cexEr2)
      text(x = (usr[1] + usr[2])/2, y=usr[3], yr2, col="white", cex=.8)
      
      plot(network, col="gray25", bg="black")
      plot(network3, col=fix3.colors, bg="black", fg="white", add=TRUE)
      points(error_pts3, pch=16, col="gray40", cex=cexEr3)
      text(x = (usr[1] + usr[2])/2, y=usr[3], yr3, col="white", cex=.8)
      
      mtext("Wenatchee 8-day mean ('C)",side =3, outer=TRUE, line=-3, col="white", cex=.9)
      mtext(paste0("Julian Day ", namey),side =1, outer=TRUE, line=-3, col="white", cex=.7)
      
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend(x=grconvertX(c(0.9,0.99), from='npc'), 
             y=grconvertY(c(0.6, 0.8), from='npc'), title="8-day Mean ('C)", fill = attr(fix3.colors, "palette"), legend = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18+"), bty = "n", cex=.7, text.col="white", xpd= TRUE);
      
      
      dev.off()
    }


setwd(paste0(mainPath, longBasin, "/", yr3, "/graphics5/"))

system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png example.mpeg')



###############################################
#
# Edited 5 may 2016 to add animation output
# Edited 25 May 2016 to add HUC display
# Edited 13 June 2016 to add multiple frame display
# Edited 15 June 2016 to add annual profiles back in 
################################################

library(rgdal)
library(RColorBrewer)
library(classInt)

  yr1 <- 2012
  yr2 <- 2013
  yr3 <- 2014

ptsPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", "Mean_models")
shpPath3 <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yr3, "/", "graphic_shapes")
shpPath2 <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yr2, "/", "graphic_shapes")
shpPath1 <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yr1, "/", "graphic_shapes")

  netname <- paste0(basin, "_base_net")
  netname1 <- paste0(basin, "_", yr1, "_8D_mn")
  netname2 <- paste0(basin, "_", yr2, "_8D_mn") 
  netname3 <- paste0(basin, "_", yr3, "_8D_mn")
  
  ptsname1 <- paste0(basin, "_Error_", yr1, "_8D_Mn")
  ptsname2 <- paste0(basin, "_Error_", yr2, "_8D_Mn")
  ptsname3 <- paste0(basin, "_Error_", yr3, "_8D_Mn")

setwd(ptsPath)

  error_pts3 <- readOGR(dsn=".", layer = ptsname3)

setwd(shpPath1)

  error_pts1 <- readOGR(dsn=".", layer = ptsname1)

setwd(shpPath2)

  network2 <- readOGR(dsn=".", layer = netname2)
  network2@data <- network2@data[,-4]
  error_pts2 <- readOGR(dsn=".", layer = ptsname2)

setwd(shpPath3)
  network <- readOGR(dsn=".", layer = netname)

  network1 <- readOGR(dsn=".", layer = netname1)
  network1@data <- network1@data[,-4]

  network3 <- readOGR(dsn=".", layer = netname3)


#HUCs <- readOGR(dsn=".", "Wen_HUC5")

  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)

##################### basin means ####################
  
  means1 <- colMeans(network1@data[4:49])
  SDs1 <- colStdevs(network1@data[4:49])
  yplus1 <- means1 + SDs1
  yminus1 <- means1 - SDs1
  df1 <- data.frame(means=means1, SDs=SDs1, names=namesnum)
  
  
  means2 <- colMeans(network2@data[4:49])
  SDs2 <- colStdevs(network2@data[4:49])
  yplus2 <- means2 + SDs2
  yminus2 <- means2 - SDs2
  df2 <- data.frame(means=means2, SDs=SDs2, names=namesnum)
  
  means3 <- colMeans(network3@data[4:49])
  SDs3 <- colStdevs(network3@data[4:49])
  yplus3 <- means3 + SDs3
  yminus3 <- means3 - SDs3
  df3 <- data.frame(means=means3, SDs=SDs3, names=namesnum)
  
  fix1a <- classIntervals(means1, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix1a.colors <- findColours(fix1a,pal=seis)
  
  fix2a <- classIntervals(means2, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix2a.colors <- findColours(fix2a,pal=seis)
  
  fix3a <- classIntervals(means3, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix3a.colors <- findColours(fix3a,pal=seis)
  
  

##################### graphics out #################################

sequ <- c(1:46)
namer <- sprintf('%03d', sequ)


names.out <- colnames(network3@data[4:49])
namesnum <- as.numeric(gsub("Tmn_14_", "", colnames(network3@data[4:49])))

op <- par(mfrow = c(1,3),
          mar = c(0,0,1,1) + 0.1)
usr <- par("usr")

  for (i in 4:49)
    {
      
      namey <- gsub("Tmn_14_", "", colnames(network3@data)[i])
      
      filename <- paste0(mainPath, longBasin, "/", yr1, "/graphics3/", namer[i-3], ".png", sep="") 
      png(filename=filename, res = 300, width = 1800, height = 1500, units = "px", bg="black")
      
      fix1 <- classIntervals(network1@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix1.colors <- findColours(fix1,pal=seis)
      
      fix2 <- classIntervals(network2@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix2.colors <- findColours(fix2,pal=seis)
      
      fix3 <- classIntervals(network3@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix3.colors <- findColours(fix3,pal=seis)
      
      cexEr1 <- ifelse(abs(error_pts1@data[,i-1]) <= 1, 0.5,
                      ifelse(abs(error_pts1@data[,i-1])>1, 0.75,
                             ifelse(abs(error_pts1@data[,i-1])>2, 1.0,
                                    ifelse(abs(error_pts1@data[,i-1])>3, 1.25, NA))))

      cexEr2 <- ifelse(abs(error_pts2@data[,i-1]) <= 1, 0.5,
                       ifelse(abs(error_pts2@data[,i-1])>1, 0.75,
                              ifelse(abs(error_pts2@data[,i-1])>2, 1.0,
                                     ifelse(abs(error_pts2@data[,i-1])>3, 1.25, NA))))
      
      cexEr3 <- ifelse(abs(error_pts3@data[,i-1]) <= 1, 0.5,
                       ifelse(abs(error_pts3@data[,i-1])>1, 0.75,
                              ifelse(abs(error_pts3@data[,i-1])>2, 1.0,
                                     ifelse(abs(error_pts3@data[,i-1])>3, 1.25, NA))))
      s <- i-4
      plx1 <- function() {plot(namesnum[1:(i-3)], means1[1:(i-3)], type="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"); 
        segments(x0=namesnum[1:s], y0=means1[1:s], x1=namesnum[2:(s+1)], y1=means1[2:(s+1)], col=fix1a.colors)}
      
      plx2 <- function() {plot(namesnum[1:(i-3)], means2[1:(i-3)], type="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"); 
        segments(x0=namesnum[1:s], y0=means2[1:s], x1=namesnum[2:(s+1)], y1=means2[2:(s+1)], col=fix2a.colors)}
      
      plx3 <- function() {plot(namesnum[1:(i-3)], means3[1:(i-3)], type="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"); 
        segments(x0=namesnum[1:s], y0=means3[1:s], x1=namesnum[2:(s+1)], y1=means3[2:(s+1)], col=fix3a.colors)}
      
      
      par(mfrow = c(1,3),
                oma = c(0,0,0,0),
                mar = c(0,0,0,0),
                mgp = c(2,0,0),
                xpd=NA)
      
      
      plot(network, col="gray25", bg="black")
      plot(network1, col=fix1.colors, bg="black", fg="white", add=TRUE)
      points(error_pts1b, pch=16, col="grey92", cex=cexEr1)
      text(x = (usr[1] + usr[2])/2, y=usr[3]*1.04, yr1, col="white", cex=.8)
     
      tmp2 <- subplot(plx1(),
                      x=grconvertX(c(0.08,0.45), from='npc'), 
                      y=grconvertY(c(0.19, 0.29), from='npc'),
                      size=c(1,1.5), vadj=0.5, hadj=0.5, 
                      pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
                    
      plot(network, col="gray25", bg="black")
      plot(network2, col=fix2.colors, bg="black", fg="white", add=TRUE)
      points(error_pts2b, pch=16, col="grey92", cex=cexEr2)
      text(x = (usr[1] + usr[2])/2, y=usr[3]*1.04, yr2, col="white", cex=.8)
      
      tmp2 <- subplot(plx2(),
                      x=grconvertX(c(0.08,0.45), from='npc'), 
                      y=grconvertY(c(0.19, 0.29), from='npc'),
                      size=c(1,1.5), vadj=0.5, hadj=0.5, 
                      pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
      
      plot(network, col="gray25", bg="black")
      plot(network3, col=fix3.colors, bg="black", fg="white", add=TRUE)
      points(error_pts3, pch=16, col="grey92", cex=cexEr3)
      text(x = (usr[1] + usr[2])/2, y=usr[3]*1.04, yr3, col="white", cex=.8)
      
      tmp2 <- subplot(plx3(),
                      x=grconvertX(c(0.08,0.45), from='npc'), 
                      y=grconvertY(c(0.19, 0.29), from='npc'),
                      size=c(1,1.5), vadj=0.5, hadj=0.5, 
                      pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
      
      mtext("Wenatchee 8-day mean (°C)",side =3, outer=TRUE, line=-4, col="white", cex=.9)
      mtext("Steelhead extent",side =3, outer=TRUE, line=-6, col="white", cex=.7)
      mtext(paste0("Julian Day ", namey),side =1, outer=TRUE, line=-3, col="white", cex=.7)
      
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend(x=grconvertX(c(0.91,0.99), from='npc'), 
             y=grconvertY(c(0.6, 0.8), from='npc'), title="8-day Mean (°C)", fill = attr(fix3.colors, "palette"), legend = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18+"), bty = "n", cex=.65, text.col="white", xpd= NA);
      legend(x=grconvertX(c(0.92,0.99), from='npc'), 
             y=grconvertY(c(0.1, 0.18), from='npc'), title="Model error (°C)", legend = c("0-1","1-2","2-3","3+"), bty = "n", pch=16, pt.cex=c(0.5, 0.75, 1.0, 1.25), col = "grey92", cex=.6, text.col="white", xpd=NA);
      
      
      dev.off()
    }


setwd(paste0(mainPath, longBasin, "/", yr1, "/graphics3/"))

system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png example1.mpeg')


###########################
pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
nf <- layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), respect = TRUE)
layout.show(nf)

