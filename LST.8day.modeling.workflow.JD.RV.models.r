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

  basin <- "JD"
  midBasin <- "JohnDay"
  longBasin <- "JohnDay"
  yrPath <- 13
  yearPath <- "2013"
  dataPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/climate_Analysis/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"

##############################################################################################
# This section reads in the 1km LST data for a year, uses a 4th order polynomial to fill in Julian day 1 & 365 (if they are missing)
# then fills any remaining gaps across the year at each pixel with a linear interpolation.
# #############################################################################################
  
  
  dataDir <- "D:/OneDrive/work/GIS/8_day_1km_LST/LST_s1_2013/"
  setwd(paste0(dataDir))
  
  
  LST.in <- read.dbf(paste0("JD_13_LST.dbf"))
  GrPolID <- LST.in[,1]
  LST.in <- LST.in[,2:47]
  LST.in[LST.in<0.1] = NA
  
  LST.names<-colnames(LST.in)
  
  tLST.in <- t(LST.in)
  tLST.out <- na.spline(tLST.in)
  CLST <- tLST.out*0.02-273.15
  tCLST <- t(CLST)
  CLST.out <- as.data.frame(tCLST)
  
  CLST.out$GRID_CODE <- GrPolID
  
  colnames(CLST.out)[1:46] <- LST.names
  
  setwd(paste0(mainPath, midBasin))
  
  write.dbf(CLST.out, file = paste0("LST", yrPath, "_", basin, "_interp.dbf"))
            
# ###################################################################
#This part calculates the zonal mean of the daily LST for each RCA
# ###################################################################
            
    out_Preds <- read.dbf(paste0("LST", yrPath, "_", basin, "_interp.dbf")) 
            
    weights <- read.csv(paste0(basin, "_rca_area_wgts.csv"))
    weights[weights$area_wgts<0.01, "area_wgt"] = 0.01
            
    rcas <- unique(unlist(weights$RCA_ID))
    rca_zonal <- matrix(ncol=46, nrow=length(rcas))
            
            l <- 1
            for(i in rcas)  
            {
              pixels <- weights[weights$RCA_ID == i, "GRIDCODE"]
              wgts <- weights[weights$RCA_ID == i, "area_wgt"]
              
              for (j in 1:46)
              {
                daily <- out_Preds[out_Preds$GRID_CODE %in% pixels, j]
                zonal_mn <- weighted.mean(daily, wgts, na.rm = TRUE)
                rca_zonal[l,j] <- zonal_mn
              }
              l <- l+1
            }
            
            names.out <- sprintf("%03d", newnamesnum)
            colnames(rca_zonal)[1:46] <- names.out
            rca_zonal <- as.data.frame(rca_zonal)
            rca_zonal$RCAID <- rcas
            
            write.dbf(rca_zonal, file = paste0("LST", yrPath, "_", basin, "_RCA.dbf"))
            write.csv(rca_zonal, file = paste0("LST", yrPath, "_", basin, "_RCA.csv"))
            
            
            #    
#################################################################
#Logger data processing part
  
#################################################################
  Alldata <- data.frame(mup=NULL)
  
  setwd(paste0(mainPath, longBasin))
  
  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))
      #######OR#######
  LST.in <- read.csv(paste0("LST", yrPath, "_", basin, "_RCA.csv"))
                     
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  
  ID.in <- read.csv(paste0("All_", basin, "_CHaMP_sites_RCA_elev.csv"), stringsAsFactors=FALSE)
  colnames(ID.in)[4] <- "RCAID"
  
  setwd(paste0(mainPath, longBasin))
  
  Log.in <- read.csv(paste0(basin, "_", yearPath, "_summer_daily_logger_data.csv"), stringsAsFactors=FALSE)
  
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
        RCAID <- ID.in$RCAID[ID.in$SiteName == i]
        Elev <- ID.in$Elev_M[ID.in$SiteName == i]
        Max <- max(Log.site$DailyMax)
        Mn <- mean(Log.site$DailyMn)
        MnRng <- mean(Log.site$DailyRange)
        AbRng <- Max - min(Log.site$DailyMin)
        LST.site <- unlist(LST.in[LST.in$RCAID == RCAID,23:33])
        mnLST <- mean(LST.site)
        mxLST <- max(LST.site)
        rngLST <- max(LST.site) - min(LST.site)
        precip <- y01_15_ppt[yrPath]
        
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
  RVID <- unique(RV.in$RCA_ID)
  
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
  yRng <- Data.means$MnRng
  xRng <- Data.means$mnLST
  yMax <- Data.means$Max
  xMax <- Data.means$mxLST
  
  plot(evt, y)

  mod <- lm(y ~ x + e + evt + depr + bps + year)
  sum_mod <- summary(mod)
  
  mod2 <- lm(yMax ~ xMax + e + evt + depr + bps + year)
  sum_mod2 <- summary(mod2)
  sum_mod2
 
  mod3 <- lm(yRng ~ xRng + e + evt + depr + bps + year)
  sum_mod3 <- summary(mod3)
  sum_mod3
  
  xRng <- Data.means$rngLST
  yRng <- Data.means$AbRng
  ###################################
  
  
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

