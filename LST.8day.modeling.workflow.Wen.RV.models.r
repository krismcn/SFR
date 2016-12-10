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
library(foreign)

  basin <- "Wen"
  midBasin <- "Wenatchee"
  longBasin <- "Wenatchee"
  dataPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/climate_Analysis/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
    
#################################################################
#Logger data processing part
  
#################################################################
  Alldata <- data.frame(mup=NULL)
  
  setwd(paste0(mainPath, longBasin))
  
  ID.in <- read.csv(paste0(basin, "_sites_elev.csv"), stringsAsFactors=FALSE)
  colnames(ID.in)[5] <- "RCAID"
  colnames(ID.in)[4] <- "Elev"
  
  yrPath <- "14"
  year <- "2014"
  
  LST.in <- read.csv(paste0("LST", yrPath, "_", basin, "_RCA.csv"), stringsAsFactors = F)
                            ########### OR ###########
  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))
  
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  
  
  Log.in <- read.csv(paste0(basin, "_", year, "_summer_daily_logger.csv"), stringsAsFactors=FALSE)
  
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)
  model.data <- data.frame(mup=NULL)
  LST.Log.out <- data.frame (mup = NULL)
  
  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteName  == i,]
      Log.site <- as.data.frame(Log.site)
      
      if (dim(Log.site)[1] > 90)
      {
        RCAID <- ID.in$RCAID[ID.in$SiteID == i]
        Elev <- ID.in$Elev[ID.in$SiteID == i]
        Max <- max(Log.site$DailyMax)
        Mn <- mean(Log.site$DailyMn)
        MnRng <- mean(Log.site$DailyRange)
        AbRng <- Max - min(Log.site$DailyMin)
        LST.site <- unlist(LST.in[LST.in$RCAID == RCAID,23:33])
        mnLST <- mean(LST.site)
        mxLST <- max(LST.site)
        rngLST <- max(LST.site) - min(LST.site)
        precip <- y01_15_ppt[14]
        
        data <- data.frame(RCAID=RCAID, Elev=Elev, Max=Max, Mn=Mn, MnRng=MnRng, AbRng=AbRng, mnLST=mnLST, mxLST=mxLST, rngLST=rngLST, precip=precip, year=year)
        model.data <- rbind(model.data, data)
       
      } 
    }
  
  
  Alldata <- rbind(Alldata, model.data)
  
# ################################################################
# This section reads in the formatted data from the temperature models
# and merges it with the RV data by site.
# ################################################################
  
  setwd(paste0(mainPath, longBasin))
  ID.in <- read.csv(paste0(basin, "_sites_elev.csv"), stringsAsFactors=FALSE)
  
  setwd(dataPath)
  RV.in <- read.dbf(paste0(basin, "_RV_RCA.dbf"))

  setwd(paste0(mainPath, longBasin, "/", year, "/"))
 
  data.in <- Alldata
  #data.in <- read.csv(paste0(basin, "_", year, "_8Day_model_data.csv"), header=TRUE, stringsAsFactors=FALSE)
  

  RVID <- unique(RV.in$RV_reachID)
  LST.RV.out <- data.frame (mup = NULL)
  rcas <- unique(data.in$RCAID)
  RVID <- unique(RV.in$rca_id)
  
    for (i in RVID) 
      { 
        
        RCAID <- i
        
        if(RCAID %in% rcas)
          
          {
              data <- data.in[data.in$RCAID == RCAID,]
              data$EVT <- RV.in$EVT_MEAN[RV.in$RV_reachID == i]
              data$BPS <- RV.in$BPS_MEAN[RV.in$RV_reachID == i]
              data$DepR <- RV.in$DEP_RATIO[RV.in$RV_reachID == i]
              data$ConvC <- RV.in$CONV_CODE[RV.in$RV_reachID == i]
              data$ConvT <- RV.in$CONV_TYPE[RV.in$RV_reachID == i]
              
              LST.RV.out <- rbind(LST.RV.out, data)
          }
          
      }
  
 #  setwd(paste0(mainPath, longBasin, "/", year, "/"))
 #  write.csv(x=LST.RV.out, file=paste0(basin, "_", year, "_RV_model_data.csv"), row.names=FALSE)
 #  Data.means <- data.frame(siteName=sites)
 #  
 # for (i in sites)
 #  {
 #    siteData <- LST.RV.out[LST.RV.out$z < 263 & LST.RV.out$z > 173 & LST.RV.out$SiteName == i,]
 #    Data.means$xMn[Data.means$siteName == i] <- mean(siteData$x)
 #    Data.means$yMn[Data.means$siteName == i] <- mean(siteData$y)
 #    Data.means$elev[Data.means$siteName == i] <- siteData$e[1]
 #    Data.means$EVT[Data.means$siteName == i] <- siteData$EVT[1]
 #    Data.means$DepR[Data.means$siteName == i] <- siteData$DepR[1]
 #    Data.means$BPS[Data.means$siteName == i] <- siteData$BPS[1]
 #    Data.means$ConvC[Data.means$siteName == i] <- siteData$ConvC[1]
 #  }  
    
 Data.means <- LST.RV.out   
 write.csv(Data.means, file=paste0(basin, "RV_model_data_12_14.csv"))
# ################
# Adding in the NBCD data
# ################
 
setwd(dataPath)
 
 Data.means <- read.csv(paste0(basin, "RV_model_data_12_14.csv"))
 View(Data.means)
 setwd(paste0(dataPath, "NBCD/"))
 baw <- read.csv(paste0(basin, "_RCA_BAW.csv"))
 
 Data <- merge(Data.means, baw, by="RCAID", all.x=TRUE, all.y=FALSE)
 Data.means <- Data
# ###############################
# full year
# ##################################

  y <- Data.means$Mn
  x <- Data.means$mnLST
  e <- Data.means$Elev
  evt <- Data.means$EVT
  depr <- Data.means$DepR
  bps <- Data.means$BPS
  year <- Data.means$precip
  yMnRng <- Data.means$MnRng
  yAbsRng <- Data.means$AbRng
  xRng <- Data.means$rngLST
  yMax <- Data.means$Max
  xMax <- Data.means$mxLST
  baw <- Data.means$NBCD_BAW
  plot(evt, y)
  plot(baw, y)

  mod <- lm(y ~ x + e + evt + depr + bps + baw + year)
  sum_mod <- summary(mod)
  sum_mod
  
  mod2 <- lm(yMax ~ xMax + e + evt + depr + bps + baw + year)
  sum_mod2 <- summary(mod2)
  sum_mod2
  
  mod3 <- lm(y ~ x + evt + year)
  sum_mod3 <- summary(mod3)
  sum_mod3
  
  mod4 <- lm(y ~ x + depr + year)
  sum_mod4 <- summary(mod4)
  sum_mod4
  
  mod5 <- lm(y ~ x + bps + year)
  sum_mod5 <- summary(mod5)
  sum_mod5
  
  mod6 <- lm(y ~ x + baw + year)
  sum_mod6 <- summary(mod6)
  sum_mod6
  
  mod7 <- lm(yMnRng ~ xRng + e + evt + depr + bps + baw + year)
  sum_mod7 <- summary(mod7)
  sum_mod7
  
  mod8 <- lm(yMnRng ~ xRng + evt + year)
  sum_mod8 <- summary(mod8)
  sum_mod8
  
  mod9 <- lm(yMnRng ~ xRng + depr + year)
  sum_mod9 <- summary(mod9)
  sum_mod9
  
  mod10 <- lm(yMnRng ~ xRng + bps + year)
  sum_mod10 <- summary(mod10)
  sum_mod10
  
  mod11 <- lm(yMnRng ~ xRng + baw + year)
  sum_mod11 <- summary(mod11)
  sum_mod11
  
  mod12 <- lm(yMax ~ xMax + evt + year)
  sum_mod12 <- summary(mod12)
  sum_mod12
  
  mod13 <- lm(yMax ~ xMax + depr + year)
  sum_mod13 <- summary(mod13)
  sum_mod13
  
  mod14 <- lm(yMax ~ xMax + bps + year)
  sum_mod14 <- summary(mod14)
  sum_mod14
  
  mod15 <- lm(yMax ~ xMax + baw + year)
  sum_mod15 <- summary(mod15)
  sum_mod15
  
  mod16 <- lm(yAbsRng ~ xRng + e + evt + depr + bps + baw + year)
  sum_mod16 <- summary(mod16)
  sum_mod16
  
  mod17 <- lm(yAbsRng ~ xRng + evt + year)
  sum_mod17 <- summary(mod17)
  sum_mod17
  
  mod18 <- lm(yAbsRng ~ xRng + depr + year)
  sum_mod18 <- summary(mod18)
  sum_mod18
  
  mod19 <- lm(yAbsRng ~ xRng + bps + year)
  sum_mod19 <- summary(mod19)
  sum_mod19
  
  mod20 <- lm(yAbsRng ~ xRng + baw + year)
  sum_mod20 <- summary(mod20)
  sum_mod20
  
  
  
  pred.y <- predict(mod)
  plot(pred.y, y, main = "8-day Mean Full Year")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  gvmodel <- gvlma(mod)
  summary(gvmodel)
  plot(mod, which= 1:6)
  outlierTest(mod)
  qqPlot(mod, main="QQ Plot Full Year")
  spreadLevelPlot(mod)
  plot(pred.y, mod$residuals, main="Model diagnostics Full Year", xlab="Predicted", ylab="Residuals")
  
  pressstat_sum <- PRESS(mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
  
  library(pls)
  mod2 <- plsr(y ~ x + I(x^2) + e, validation = "LOO")
  p2 <- R2(mod2)
  detach("package:pls", unload=TRUE)
  
  
  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "full year"
  pred.out[,5] <- yearPath
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_Mean_", basin, "_", yearPath, "_full_year.csv"),sep = ",", col.names=F)  
  
  plot(pred.out[,1], pred.out[,2])
  summer_pred <- subset(pred.out, z > 181 & z < 258)
  points(summer_pred[, 1], summer_pred[,2], pch = 16, col = "green")
  abline(0,1)

  fit <- lm(y~ pred.y)
  plot(fit)
  summary(fit)
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

############################################
# These parts are used as needed
############################################

  new_mod <- lm(y ~ I(x^2) + z + e)
  new_sum <- summary(new_mod)
  new_sum
  pressstat_new <- PRESS(new_mod, verbose = "FALSE")
  pressstat_new$stat

  mod <- new_mod
###########################################
    
  
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  
  pred.y[pred.y < -0.5] = -0.5
  
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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_Mn_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  
  
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

############################################
# These parts are used as needed
############################################

  new_mod <- lm(y ~ x + z + e)
  new_sum <- summary(new_mod)
  new_sum
  pressstat_new <- PRESS(new_mod, verbose = "FALSE")
  pressstat_new$stat
  
  mod <- new_mod
###########################################

  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5
  plot(pred.y, y, main = "8-day Mn Fall Leg")
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
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_Mn_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  

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

  pred.y <- read.csv(paste0("jk_pred_v_y_Mn_", basin, "_", yearPath, "_sp_fall.csv"), stringsAsFactors = FALSE)
  colnames(pred.y) <- c("Y", "PredY", "JulDay", "Season", "Year")
  
  plot(pred.y$PredY, pred.y$Y, pch=16, col="blue", main=paste0("Mn 8-day stream temp ", basin, " ", yearPath), xlab="Predicted", ylab="Observed")
  abline(0,1)
  abline(lm(pred.y$Y~ pred.y$PredY), col="blue")
  fit <- lm(pred.y$Y~ pred.y$PredY)
  plot(fit)
  summary(fit)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day temp estimates for 1 July - 30 Sept
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

  NoNA.xyz <- read.csv(paste0(basin, "_", yearPath, "_8Day_model_data.csv"), stringsAsFactors = TRUE)
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


  modelPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", yearPath, "_temp_CHaMP/")
  ptsPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", "Mean_models")
  shpPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", midBasin, "/", yearPath, "/", "graphic_shapes")

  netname <- paste0(basin, "_", yearPath, "_8D_mn")
  netname2 <-paste0(basin, "_2013_8D_mn") 
  ptsname <- paste0(basin, "_Error_", yearPath, "_8D_mn")

setwd(ptsPath)

  error_pts <- readOGR(dsn=".", "Wen_Error_2014_8D_Mn")

setwd(shpPath) 

  network <- readOGR(dsn=".", layer = netname)

  network2 <- readOGR(dsn=".", layer = netname2)
  network2@data <- network2@data[,-4]
  
  HUCs <- readOGR(dsn=".", "Wen_HUC5")
  
  means2 <- colMeans(network2@data[4:49])
  SDs2 <- colStdevs(network2@data[4:49])
  yplus2 <- means2 + SDs2
  yminus2 <- means2 - SDs2
  df2 <- data.frame(means=means, SDs=SDs2, names=namesnum)
  sequ <- c(1:46)
  namer <- sprintf('%03d', sequ)
  fix4 <- classIntervals(means, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix4.colors <- findColours(fix4,pal=seis)
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)
  
  names.out <- colnames(network@data[4:49])
  namesnum <- as.numeric(gsub("Tmn_14_", "", colnames(network@data[4:49])))

  op <- par(mfrow = c(1,3),
            mar = c(0,0,1,1) + 0.1)

  for (i in 4:49)
    {
      
      namey <- gsub("Tmn_14_", "", colnames(network@data)[i])
      
      filename <- paste0(mainPath, longBasin, "/", yearPath, "/graphics4/", namer[i-3], ".png", sep="")
      png(filename=filename, res = 300, width = 1500, height = 1500, units = "px", bg="black")
      
      fix3 <- classIntervals(network@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
      fix3.colors <- findColours(fix3,pal=seis)
      
      op <- par(mfrow = c(1,2),
                oma = c(0,2,0,0),
                mar = c(1,1,1,1.5) + 0.1,
                mgp = c(2,0,0),
                xpd=NA)
      plot(network, col=fix3.colors, bg="black", fg="white")
      text(x = (usr[1] + usr[2])/2, y=usr[3], yearPath, col="white", cex=.8)
      
      plot(network2, col=fix3.colors, bg="black", fg="white")
      text(x = (usr[1] + usr[2])/2, y=usr[3], "2013", col="white", cex=.8)
      
      legend("right", title="8-day Mean ('C)", fill = attr(fix3.colors, "palette"), legend = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18+"), bty = "n", cex=.5, inset=c(-.2,0), text.col="white");
      mtext("Wenatchee 8-day mean ('C)",side =3, outer=TRUE, line=-3, col="white", cex=.9)
      mtext(paste0("Julian Day ", namey),side =1, outer=TRUE, line=-3, col="white", cex=.7)
      
      
      dev.off()
    }


setwd(paste0(mainPath, longBasin, "/", yearPath, "/graphics4/"))

system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png example.mpeg')


###########################
pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
nf <- layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), respect = TRUE)
layout.show(nf)

