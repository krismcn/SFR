############################################################################################################
# This set of R scripts processes logger data for predictive stream temperature modeling
# Created: 21 Feb 2017
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

  basin <- "USal"
  midBasin <- "USalmon"
  longBasin <- "UpperSalmon"
  yrPath1 <- "14"
  yrPath2 <- "15"
  dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  epPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/"

  setwd(paste0(mainPath, basin))
  
########### NorWeST logger data processing ###########
  
  setwd(paste0(epPath, midBasin))
  
  norWest <- read.csv("NorWeST_logger.csv", stringsAsFactors = FALSE)
  
  setwd(paste0(mainPath, basin))
  
  full_year <- read.csv("NorWeST_full_year_logger_data.csv", stringsAsFactors = FALSE)
  
  cbs_sites <- read.csv("USal_NorWest_CBS_sites.csv", stringsAsFactors = FALSE)
  
  full_year_sites <- read.csv("USal_NorWest_full_year_sites.csv", stringsAsFactors = FALSE)
  
  
  full_yr_logger <- full_year[full_year$OrigID %in% full_year_sites$Orig_ID, ]
  
  full_yr_logger <- merge(full_yr_logger, full_year_sites, by.x = "OrigID", by.y="Orig_ID")
  
  USal_NorWesT <- norWest[norWest$OBSPRED_ID %in% cbs_sites$OBSPRED_ID,]
  
  USal_NorWesT <- merge(USal_NorWesT, cbs_sites, by.x="OBSPRED_ID", by.y="OBSPRED_ID")
  
  USal <- full_yr_logger[full_yr_logger$subBasin == "USal",]
  
  Pahs <- full_yr_logger[full_yr_logger$subBasin == "Pahsimeroi",]
  
  Panther <- full_yr_logger[full_yr_logger$subBasin == "MSal_Panther",]
  
  Lemhi <- full_yr_logger[full_yr_logger$subBasin == "Lemhi",]
  
 
  USal_summer <- USal_NorWesT[USal_NorWesT$subBasin == "USal",]
  
  Pahs_summer <- USal_NorWesT[USal_NorWesT$subBasin == "Pahsimeroi",]
  
  Panther_summer <- USal_NorWesT[USal_NorWesT$subBasin == "MSal_Panther",]
  
  Lemhi_summer <- USal_NorWesT[USal_NorWesT$subBasin == "Lemhi",]
  
  plot(USal$JulDay, USal$DailyMean, main="USal")
  unique(USal$OrigID)
  points(USal_summer$JulDay, USal_summer$DailyMean, pch=16, col="blue")
  
  plot(Pahs$JulDay, Pahs$DailyMean, main="Pahsimeroi")
  unique(Pahs$OrigID)
  unique(Pahs$Year)
  points(Pahs_summer$JulDay, Pahs_summer$DailyMean, pch=16, col="blue")
  
  plot(Panther$JulDay, Panther$DailyMean, main = "Panther")
  unique(Panther$OrigID)
  unique(Panther$Year)
  points(Panther_summer$JulDay, Panther_summer$DailyMean, pch=16, col="blue")
  
  plot(Lemhi$JulDay, Lemhi$DailyMean, main = "Lemhi")
  unique(Lemhi$OrigID)
  unique(Lemhi$Year)
  points(Lemhi_summer$JulDay, Lemhi_summer$DailyMean, pch=16, col="blue")
  setwd(paste0(mainPath, "USal/Lemhi/"))
  write.csv(Lemhi, "Lemhi_CBS_data.csv")
  write.csv(Lemhi_summer, "Lemhi_NorWeST_data.csv")
  
  
#######################################################################
# This part reads the LST data and gap-fills using a cubic spline function
#######################################################################
  
  basin <- "Lem"
  midBasin <- "Lemhi"
  longBasin <- "Lemhi"
  yrPath <- "11"
  yearPath <- 2011
  dataPath <- "D:/OneDrive/work/GIS/8_day_1km_LST/LST_s1_"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  
  
  setwd(paste0(dataPath, yearPath, "/"))
  
  LST.in <- read.dbf(paste0(basin, "_", yrPath, "_LST.dbf"))
  PointID <- LST.in[,1]
  LST.in <- LST.in[,2:47]
  LST.in[LST.in<0.1] = NA
  
  LST.names<-colnames(LST.in)
  LST.names <- order(LST.names)
  LST.in<- LST.in[,(LST.names)]
  
  tLST.in <- t(LST.in)
  tLST.out <- interpNA(tLST.in, method = "linear")
  CLST <- tLST.out*0.02-273.15
  tCLST <- t(CLST)
  CLST.out <- as.data.frame(tCLST)
  
  
  myFunc <- function(x)
    {
      if (is.na(x[1]))
      {
        x[1] <- x[2]
      }
      x
    }
  
  LST.filled <- apply(CLST.out, 1, function(x) myFunc(x))
  
  myFunc2 <- function(x)
    {
      if (is.na(x[45]))
      {
        x[45] <- x[44]
      }
      x
    }
  
  tLST.filled <- t(LST.filled)
  LST.filled <- apply(tLST.filled, 1, function(x) myFunc2(x))
  
  myFunc3 <- function(x)
    {
      if (is.na(x[46]))
      {
        x[46] <- x[45]
      }
      x
    }
  
  tLST.filled <- t(LST.filled)
  LST.filled <- apply(tLST.filled, 1, function(x) myFunc3(x))
  tLST.filled <- t(LST.filled)
  LST.out <- as.data.frame(tLST.filled)
  LST.out$PointID <- PointID
  
  plot(1:46, LST.out[1,1:46])
  plot(1:46, LST.out[1000,1:46])
  
  setwd(paste0(mainPath, longBasin, "/"))
  write.dbf(LST.out, file = paste0("LST", yrPath, "_", basin, "_interp.dbf"))
  
  
####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################
  
  setwd(paste0(mainPath, longBasin, "/"))
  out_Preds <- read.dbf(paste0("LST", yrPath, "_", basin, "_interp.dbf"))
  
  weights <- read.dbf(paste0(basin, "_rca_area_wgts.dbf"))
  weights[weights$area_wgts<0.01, "area_wgts"] = 0.01
  
  rcas <- unique(unlist(weights$RCAID))
  rca_zonal <- matrix(ncol=46, nrow=length(rcas))
  
  
  l <- 1
  
  for(i in rcas)	
    {
      pixels <- weights[weights$RCAID == i, "GRIDCODE"]
      wgts <- weights[weights$RCAID == i, "area_wgts"]
      
      for (j in 1:46)
      {
        daily <- out_Preds[out_Preds$PointID %in% pixels, j]
        zonal_mn <- weighted.mean(daily, wgts, na.rm = TRUE)
        rca_zonal[l,j] <- zonal_mn
      }
      l <- l+1
    }
  
  
  colnames(rca_zonal)[1:46] <- colnames(out_Preds)[1:46]
  rca_zonal <- as.data.frame(rca_zonal)
  rca_zonal$RCAID <- rcas
  
  plot(1:46, rca_zonal[1,1:46])
  plot(1:46, rca_zonal[1000,1:46])
  write.dbf(rca_zonal, file = paste0("LST", yrPath, "_", basin, "_RCA.dbf"))
  
  
  
#################################################################
#Logger prediction modeling part
#First summarizes the daily data into 8-day means
#################################################################
  
  setwd(paste0(mainPath, "USal/", longBasin))
 
  Log.in <- read.csv(paste0(basin, "_CBS_NorWeST_", yearPath, "_logger_data.csv"), stringsAsFactors = FALSE)
  
  
 
  dates <- matrix(nrow=365, ncol=1)
  dates[,1] <- 1:365
  colnames(dates) <- c("Days")
  
  SiteID <- unique(Log.in$SiteID)
  SiteID <- as.matrix(SiteID)
  Log.8Day.out <- data.frame (mup = NULL)
  
  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteID == i,]
      full.year <- merge(dates, Log.site, by.x = "Days", by.y = "JulDay", all.x = TRUE)
      eightday <- rollapply(full.year$DailyMean, width = 8, FUN = mean, na.rm = T, fill = NA, align = c("left"))
      eightday <- as.matrix(eightday)
      full.year$Mn8D <- eightday
      full.year$Mn8D[361] <- mean(full.year$DailyMean[361:365])
      Log.8Day.out <- rbind(Log.8Day.out, full.year)
      
    }

  
  write.csv(Log.8Day.out, paste0(basin, "_", yearPath, "_CBS_NorWeST_8Day_Mn.csv"), row.names = FALSE)
  
  #################################################################
  #Logger prediction modeling part
  #First summarizes the daily data into 8-day maxes
  #################################################################
  
  var <- "Mx"
  setwd(paste0(mainPath, "USal/", longBasin))
  
  Log.in <- read.csv(paste0(basin, "_CBS_NorWeST_", yearPath, "_logger_data.csv"), stringsAsFactors = FALSE)
  
  
  
  dates <- matrix(nrow=365, ncol=1)
  dates[,1] <- 1:365
  colnames(dates) <- c("Days")
  
  SiteID <- unique(Log.in$SiteID)
  SiteID <- as.matrix(SiteID)
  Log.8Day.out <- data.frame (mup = NULL)
  
  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteID == i,]
      full.year <- merge(dates, Log.site, by.x = "Days", by.y = "JulDay", all.x = TRUE)
      eightday <- rollapply(full.year$DailyMaximum, width = 8, FUN = mean, na.rm = T, fill = NA, align = c("left"))
      eightday <- as.matrix(eightday)
      full.year$Mx8D <- eightday
      full.year$Mx8D[361] <- mean(full.year$DailyMaximum[361:365])
      Log.8Day.out <- rbind(Log.8Day.out, full.year)
      
    }
  
  ind <- apply(Log.8Day.out, 1, function(x) !any(is.na(x)))
  Log.8Day.out <- Log.8Day.out[ind,]
  
  write.csv(Log.8Day.out, paste0(basin, "_", yearPath, "_CBS_NorWeST_8Day_", var, ".csv"), row.names = FALSE)
  
  ################## Do some data set melding outside of R #######################
  
  Log.in <- read.csv(paste0(basin, "_", yearPath, "_all_sites_", var, "_data.csv"), stringsAsFactors = FALSE)
  
  ID.in <- read.csv(paste0(longBasin, "_sites_rca_elev.csv"), stringsAsFactors = FALSE)  
  
  setwd(paste0(mainPath, longBasin, "/"))

  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(colnames(LST.in)[1:46]);
  
  LST.Log.out <- data.frame (mup = NULL)
  SiteID <- unique(Log.in$SiteName)
  
  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteName  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$RCAID[ID.in$SiteName == i]
      Elev <- ID.in$Elev[ID.in$SiteName == i]
      
      LST.site <- matrix(ncol=3, nrow=46)
      LST.site[,1] <- as.numeric(unlist(colnames(LST.in)[1:46]))
      LST.site[,2] <- unlist(LST.in[RCAID,1:46])
      LST.site <- data.frame(LST.site)
      colnames(LST.site) <- c("JulDay", "LST", "Elev")
      LST.site[3] <- Elev
      LST.Log.site <- merge(LST.site, Log.site, by = "JulDay",all.x=TRUE, all.y = FALSE)
      
      LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
    }
  

  ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
  NoNA.xyz <- LST.Log.out[ind,]
  
  setwd(paste0(mainPath, "USal/", longBasin))
  
  NoNA.xyz <- NoNA.xyz[,c(7, 2, 1, 3, 4)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  
  write.csv(x=NoNA.xyz, file= paste0(basin, "_", yearPath, "_8Day_", var, "_model_data.csv"), row.names = FALSE)
  
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  
  plot(NoNA.xyz$x, NoNA.xyz$y)
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:(maxrow-1),]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]
  
  points(data.sp$x, data.sp$y, pch=16, col="cadetblue3")
  points(data.fall$x, data.fall$y, pch=16, col="chocolate")
  
  #####################################
  # spring/fall
  ######################################
  
  
  coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bJul=numeric(2), bElev=numeric(2))
  metrics_out <- data.frame(r2=numeric(2), RMSE=numeric(2), p2=numeric(2), RMSEP=numeric(2), N_Sites=numeric(2), N=numeric(2))
  rownames(metrics_out) <- c("Spring", "Fall")
  rownames(coeffs_out) <- c("Spring", "Fall")
  
  y <- data.sp$y
  x <- data.sp$x
  z <- data.sp$z
  e <- data.sp$e
  plot(z, y)  
  plot(x, y, main="spring")
  
  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  sum_mod
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  pressstat_sum$stat
  
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  
  pred.y[pred.y < -0.5] = -0.5
  
  plot(pred.y, y, main = "8-day Mn Spring Leg")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  gvmodel <- gvlma(mod)
  summary(gvmodel)
  #plot(mod, which= 1:5)
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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", var, "_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  
  
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
  
 
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5
  plot(pred.y, y, main = "8-day Mn Fall Leg")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  
  #plot(mod, which = 1:5)
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
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", var, "_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  
  
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
  
  setwd(paste0(mainPath, "USal/", longBasin, "/", yearPath))
  
  write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_", var, ".csv"), sep = ",", col.names=T)
  
  write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_", var, ".csv"), sep = ",", col.names=T)
  
  pred.y <- read.csv(paste0("jk_pred_v_y_", var, "_", basin, "_", yearPath, "_sp_fall.csv"), stringsAsFactors = FALSE)
  colnames(pred.y) <- c("Y", "PredY", "JulDay", "Season", "Year")
  
  plot(pred.y$PredY, pred.y$Y, pch=16, col="blue", main=paste0(var, " 8-day stream temp ", basin, " ", yearPath), xlab="Predicted", ylab="Observed")
  abline(0,1)
  abline(lm(pred.y$Y~ pred.y$PredY), col="blue")
  fit <- lm(pred.y$Y~ pred.y$PredY)
  
  summary(fit)
  
  ########################################################################################################
  # This part applies the model coefficients to the LST to generate 8-day temp estimates 
  ########################################################################################################
  setwd(paste0(mainPath, longBasin))  
  
  elev.in <- read.csv(paste0(basin, "_RCA_elev.csv"), stringsAsFactors = F)
  
  setwd(paste0(mainPath, longBasin, "/"))
  
  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))
  
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  
  LST.elev <- merge(LST.in, elev.in, by = "RCAID")
  
  
  setwd(paste0(mainPath, "USal/", longBasin, "/", yearPath))
  
  coeffs.in <- read.csv(paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_", var, ".csv"), stringsAsFactors=FALSE)
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
  varName <- paste0("T", var, "_", yrPath)
  names.out <- sprintf("%s_%03d", varName, namesnum)
  colnames(LogPred.out)[1:46] <- names.out[1:46]
  plot(namesnum, LogPred.out[1,1:46])
  points(namesnum, LogPred.out[10, 1:46], pch=16, col="blue")
  points(namesnum, LogPred.out[100, 1:46], pch=16, col="green")
  points(namesnum, LogPred.out[1000, 1:46], pch=16, col="red")
  points(namesnum, LogPred.out[2000, 1:46], pch=16, col="lightblue")
  
  LogPred.out[LogPred.out < -0.5] = -0.5
  
  write.dbf(LogPred.out, file = paste0("predt", yearPath, "_", basin, "_8D_", var, ".dbf")) 
  
##########################
#This parts formats the error by day/site info
###########################
  
  NoNA.xyz <- read.csv(paste0(basin, "_", yearPath, "_8Day_", var, "_model_data.csv"), stringsAsFactors = TRUE)
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:(maxrow-1),]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]
  
  y <- data.sp$y
  x <- data.sp$x
  z <- data.sp$z
  e <- data.sp$e
  mod <- lm(y ~ x + I(x^2) + z + e)
  pred.new <- predict(mod, data = data.sp)
  pred.new[pred.new < -0.5] = -0.5
  data.sp$pred <- unlist(pred.new)
  
  y <- data.fall$y
  x <- data.fall$x
  z <- data.fall$z
  e <- data.fall$e
  mod <- lm(y ~ x + I(x^2) + z + e)
  pred.new <- predict(mod, data=data.fall)
  pred.new[pred.new < -0.5] = -0.5
  data.fall$pred <- unlist(pred.new)
  
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
  
  write.dbf(Error.pts.out, file = paste0("Error", yearPath, "_", basin, "_8D_", var, ".dbf")) 
  write.csv(Error.pts.out, file = paste0("Error", yearPath, "_", basin, "_8D_", var, ".csv"), row.names = FALSE)
  
###############################################
# Attach the above output to a point shapefile of basin sites
# Then run for animation output
# Edited 5 may 2016 to add animation output
 
 ################################################
  
  library(rgdal)
  library(RColorBrewer)
  library(classInt)
  
  
  modelPath <- paste0(mainPath, "USal/", longBasin, "/", yearPath)
 
  netname <- paste0(basin, "_", yearPath, "_8D_", var)
  ptsname <- paste0(basin, "_Error_", yearPath, "_8D_", var)
  
  
  setwd(modelPath)
  
  error_pts <- readOGR(dsn=".", ptsname)
  
  
  network <- readOGR(dsn=".", layer = netname)
  
  
  error_pts <- spTransform(error_pts, proj4string(network))
  error_pts@data <- error_pts@data[,-3]
  
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)
  
  names.out <- colnames(network@data[4:49])
  namesnum <- as.numeric(gsub("TMx_11_", "", colnames(network@data[4:49])))
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
      
      namey <- gsub("TMx_11_", "", colnames(network@data)[i])
      
      filename <- paste0(modelPath, "/graphics/", namer[i-3], ".png", sep="")
      png(filename=filename, res = 300, width = 1500, height = 1500, units = "px", bg="black")
      
      fix3 <- classIntervals(network@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix3.colors <- findColours(fix3,pal=seis)
      
      cexEr <- ifelse(abs(error_pts@data[,i]) <= 1, 0.5,
                      ifelse(abs(error_pts@data[,i])>1, 0.75,
                             ifelse(abs(error_pts@data[,i])>2, 1.0,
                                    ifelse(abs(error_pts@data[,i])>3, 1.25, NA))))
      
      plot(network, col=fix3.colors, bg="black", fg="white")
      points(error_pts, pch=16, col="gray40", cex=cexEr)
      
      legend("topright", fill = attr(fix3.colors, "palette"), legend = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18+"), bty = "n", cex=.5, inset=c(.05,0.1), text.col="white");
      legend("topright", pch=16, col="gray40", pt.cex=c(0.5, 0.75, 1.0, 1.25), title="Model error (°C)", legend = c("0-1","1-2","2-3","3+"), bty = "n", cex=.5, inset=c(0.05,0.5), text.col="white");
      
      
      title("Lemhi 8-day max 2011 (°C)", sub = paste0("Julian Day ", namey), line=-0.9, adj=.80, col.main="white", col.sub="white", outer=FALSE, cex.main=0.5, cex.sub=0.5)
      tmp2 <- subplot(
        plot(namesnum[1:(i-2)], means[1:(i-2)], col=fix4.colors, pch=16, bty="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean (+/-SD)", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"), 
        x=grconvertX(c(0.1,0.45), from='npc'), 
        y=grconvertY(c(0.05, 0.20), from='npc'),
        size=c(1,1.5), vadj=0.5, hadj=0.5, 
        pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
      
      
      dev.off()
    }
    
  
  setwd(paste0(modelPath, "/graphics/"))
  
  system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png Lemhi_2011_8D_mx.mpeg')
  
  