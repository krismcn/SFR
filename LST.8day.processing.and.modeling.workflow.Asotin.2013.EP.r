############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for a summer max temps MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# 
# Edited Aug 2014 to add the PRESS stastic output
# Should look at PRESS output to screen to chose model - only that model will output resids and preds.

# Edited Oct 2015 for EP watershed summer output
# Edited Jan 2016 to update the gap-filling interpolation functions
# Edited Aug 2014 to add the PRESS stastic output
# Edited 28 March 2016 for Secesh

##############################################################################################
# This section reads in the 1km LST data for a year, uses a 4th order polynomial to fill in Julian day 1 & 365 (if they are missing)
# then fills any remaining gaps across the year at each pixel with a linear interpolation.
##############################################################################################

# The Secesh is a sub-basin of the SFS, so some of the file paths point there

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

  basin <- "Aso"
  midBasin <- "Asotin"
  longBasin <- "Asotin"
  yrPath <- "13"
  yearPath <- "2013"
  subDir <- "EP_temp/"
  dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  
  
           
#################################################################
#Logger prediction modeling part

#################################################################
  
  
  setwd(paste0(mainPath, midBasin))
  
  
  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))

  colnames(LST.in)<-gsub("X13", "", colnames(LST.in))

  newnamesnum <- as.numeric(colnames(LST.in)[1:46])

  setwd(paste0(mainPath, longBasin))
        
  ID.in <- read.csv(paste0(basin, "_sites_elev_rca.csv"), stringsAsFactors=FALSE)
  
  Log.in <- read.csv(paste0(basin, "_temp_data_", yearPath, ".csv"), stringsAsFactors=FALSE)
  
  setwd(paste0(mainPath, subDir, longBasin, "/", yearPath, "/"))
  
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)
  
  LST.Log.out <- data.frame (mup = NULL)
  
    for (i in SiteID) 
      { 
        Log.site <- Log.in[Log.in$SiteName  == i,]
        Log.site <- as.data.frame(Log.site)
        
        
        RCAID <- ID.in$rca_id[ID.in$SiteName == i]
        Elev <- ID.in$Elev_M[ID.in$SiteName == i]
        
        LST.site <- matrix(ncol=3, nrow=46)
        LST.site[,1] <- newnamesnum
        LST.site[,2] <- unlist(LST.in[LST.in$RCAID == RCAID,1:46])
        LST.site <- data.frame(LST.site)
        colnames(LST.site) <- c("JulDay", "LST", "Elev")
        LST.site[3] <- Elev
        LST.Log.site <- merge(LST.site, Log.site, by.x = "JulDay", by.y = "JulDay", all.x=TRUE, all.y = FALSE)
        
        LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
      }
  
  
  ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
  NoNA.xyz <- LST.Log.out[ind,]
  
  NoNA.xyz <- NoNA.xyz[,c(12, 2, 1, 3, 5)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  
  plot(NoNA.xyz$z, NoNA.xyz$y)
  plot(NoNA.xyz$x, NoNA.xyz$y)

  write.csv(x=NoNA.xyz, file=paste0(basin, "_", yearPath, "_8Day_model_data.csv"), row.names = FALSE)
  
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:(maxrow-1),]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]

  points(data.sp$x, data.sp$y, pch=16, col="cadetblue3")
  points(data.fall$x, data.fall$y, pch=16, col="chocolate")

    


################################
# full year
###################################

  setwd(paste0(mainPath, subDir, longBasin, "/", yearPath))

  y <- NoNA.xyz$y
  x <- NoNA.xyz$x
  z <- NoNA.xyz$z
  e <- NoNA.xyz$e
  plot(x, y)
  
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:maxrow,]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]


  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  pred.y <- predict(mod)
  plot(pred.y, y, main = "8-day Max Full Year")
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
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_Max_", basin, "_", yearPath, "_full_year.csv"),sep = ",", col.names=F)  
  
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
  setwd(paste0(mainPath, subDir, longBasin, "/", yearPath))

  coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bJul=numeric(2), bElev=numeric(2))
  metrics_out <- data.frame(r2=numeric(2), RMSE=numeric(2), p2=numeric(2), RMSEP=numeric(2), N_Sites=numeric(2), N=numeric(2))
  rownames(metrics_out) <- c("Spring", "Fall")
  rownames(coeffs_out) <- c("Spring", "Fall")

# 
#   data.sp <- data.sp[-131,]
#   data.sp <- data.sp[-125,]
#   data.sp <- data.sp[-119,]
#   

  y <- data.sp$y
  x <- data.sp$x
  z <- data.sp$z
  e <- data.sp$e
  plot(z, y)  
  plot(x, y)
  
  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  
  pred.y[pred.y < -0.5] = -0.5
  
  plot(pred.y, y, main = "8-day Max Spring Leg")
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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_Max_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  
  
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
  plot(x, y)

  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5
  plot(pred.y, y, main = "8-day Min Fall Leg")
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

  fall_pred <- subset(pred.out, z > 181 & z < 258)
  origY <- fall_pred[,1]
  predY <- fall_pred[,2]
  points(fall_pred[, 1], fall_pred[,2], pch = 16, col = "purple")


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

write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mx.csv"), sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_Mx.csv"), sep = ",", col.names=T)

  pred.y <- read.csv(paste0("jk_pred_v_y_Max_", basin, "_", yearPath, "_sp_fall.csv"), stringsAsFactors = FALSE)
  colnames(pred.y) <- c("Y", "PredY", "JulDay", "Season", "Year")
  
  plot(pred.y$PredY, pred.y$Y, pch=16, col="blue", main=paste0("Max 8-day stream temp ", yearPath), xlab="Predicted", ylab="Observed")
  abline(0,1)
  abline(lm(pred.y$Y~ pred.y$PredY), col="blue")
  fit <- lm(pred.y$Y~ pred.y$PredY)
  plot(fit)
  summary(fit)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################
  
  
setwd(paste0(mainPath, longBasin, "/"))

  elev.in <- read.csv(paste0(basin, "_rca_elev.csv"), stringsAsFactors=F)
  

  LST.in <- read.dbf(paste0("LST", yrPath, "_", basin, "_RCA.dbf"))

  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])

  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCAID")
  colnames(LST.elev)[1] <- "RCAID"

setwd(paste0(mainPath, subDir, longBasin, "/", yearPath))

  coeffs.in <- read.csv(paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mx.csv"), stringsAsFactors=FALSE)
  LogPred.out <- LST.elev[,c(25:36,49,1)]
  LogPred.out[,1:12] <- 0
  rcas <- unique(elev.in$RCAID)
  LST.sum <- LST.elev[,c(25:36, 49,1)]

  for (i in 1:length(rcas))  
    {
      x <- unlist(LST.sum[i,])
      maxrow <- as.numeric(which.max(x[1:12])) #either specify or let be dynamic
      midrow <- maxrow - 1
      day <- as.numeric(colnames(LST.sum)[maxrow])
      
      j <- 185
      for (l in 1:midrow)
      {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
       j <- j + 8}
      k <- day
      for (l in maxrow:12)     
      {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[13] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
       k <- k + 8}
      if (maxrow > 3)
        {
        x[(maxrow)] <- NA
        fill <- na.spline(x[1:12],)
        x[(maxrow-3):(maxrow+3)] <- fill[(maxrow-3):(maxrow+3)]
        }
      LogPred.out[i,1:12] <- x [1:12] 
    }


  LogPred.out <- as.data.frame(LogPred.out)

  LogPred.out$Basin_RCA <- paste0(basin, "_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:12]))
  varName <- paste0("Tmx_", yrPath)
  names.out <- sprintf("%s_%03d", varName, namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]
  plot(namesnum, LogPred.out[1,1:12])

write.dbf(LogPred.out, file = paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf")) 



#________________________________________________________
# summary for the summer max EPs linear fill

############################
# 15Jul-31Aug ------------------------------------------------------------------
#############################

  setwd(paste0(mainPath, subDir, longBasin))

  Au.in <- read.csv(paste0(basin, "_RCA_AUs_IP.csv"), stringsAsFactors = FALSE)
  
  
  setwd(paste0(mainPath, subDir, longBasin, "/", yearPath))

  
  Max.in <- read.dbf(paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf"))
  colnames(Max.in)<-gsub(paste0("Tmx_", yrPath, "_"), "", colnames(Max.in))

  Max.na <- matrix(ncol=49, nrow=nrow(Max.in))
  Max.na <- as.data.frame(Max.na)  
  colnames(Max.na)[1:49] <- 201:249
  
  colnames(Max.in)[1:12]
  Max.na$"201" <- Max.in$"201"
  Max.na$"209" <- Max.in$"209"
  Max.na$"217" <- Max.in$"217"
  Max.na$"225" <- Max.in$"225"
  Max.na$"233" <- Max.in$"233"
  Max.na$"241" <- Max.in$"241"
  Max.na$"249" <- Max.in$"249"
  tMax.na <- t(Max.na)
  Max.filled <- na.spline(tMax.na)
  Max.out <- t(Max.filled)
  Max.out <- as.data.frame(Max.out)
  Max.out$RCAID <- Max.in$RCAID
  colnames(Max.out)[1:49] <- 201:249

  plot(1:49, Max.out[1, 1:49],)
  points(1:49, Max.na[1, 1:49], pch=16, col="blue")



  Max.au <- merge(Max.out, Au.in, by.x = "RCAID", by.y = "RCAID")




#_________all_______________________  
rcas <- unique(Max.in$RCAID)

SumSumm.out <- data.frame ("RCAID" = rcas)


# count of days in exceedence -----------------------------------------


  for (i in 1:length(rcas)) 
    { 
      l <- rcas[i]
      MaxRCA <- Max.out[Max.out$RCAID == l,] #grab days for one RCA 
      DaysAbove12 <- length(which(MaxRCA[1,2:44]> 12)) #finds how many day in the 20July-31Aug window exceed threshold 
      DaysAbove13 <- length(which(MaxRCA[1,2:44]> 13))
      DaysAbove16 <- length(which(MaxRCA[1,2:44]> 16))
      DaysAbove18 <- length(which(MaxRCA[1,2:44]> 18))
      DaysAbove20 <- length(which(MaxRCA[1,2:44]> 20))
      DaysAbove22 <- length(which(MaxRCA[1,2:44]> 22))
      MaxMax <- max(MaxRCA[1,2:44])
      SDMax <- sd(MaxRCA[1,2:44])
      MeanMax <- mean(unlist(MaxRCA[2:44]))
      SumSumm.out$PctDays12[i] <- DaysAbove12/43
      SumSumm.out$PctDays13[i] <- DaysAbove13/43
      SumSumm.out$PctDays16[i] <- DaysAbove16/43
      SumSumm.out$PctDays18[i] <- DaysAbove18/43
      SumSumm.out$PctDays20[i] <- DaysAbove20/43
      SumSumm.out$PctDays22[i] <- DaysAbove22/43
      SumSumm.out$MxMx[i] <- MaxMax 
      SumSumm.out$SdMn[i] <- SDMax
      SumSumm.out$MnMx[i] <- MeanMax
    } 


  SumSumm.out[,2:7] <- round(SumSumm.out[,2:7], digits = 2)
  
  colnames(SumSumm.out) <- c("RCAID", paste0("Pct12_", yearPath),paste0("Pct13_", yearPath), paste0("Pct16_", yearPath), paste0("Pct18_", yearPath), paste0("Pct20_", yearPath), paste0("Pct22_", yearPath),  paste0("MxMx_", yearPath), paste0("sdMn_", yearPath), paste0("MnMx_", yearPath))

setwd(paste0(mainPath, subDir, longBasin, "/", yearPath))

write.dbf(SumSumm.out, file = paste0(basin,"_", yearPath, "_21Jul_31Aug_max_summary_All.dbf"))
write.csv(SumSumm.out, file = paste0(basin,"_", yearPath, "_21Jul_31Aug_max_summary_All.csv"), row.names = F)


#________________________________
# summmary of IP reaches only by AU
#______________________________________


SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "RCAID")
colnames(SumSummAU)[1:10] <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx")


######
#________steelhead____________________  
######

  popn <- "STHD"
  varName <- paste0(basin, "_", popn)
  ST_rcas <- unique(na.omit(Au.in[Au.in$STHD_IP == "Steelhead", "RCAID"]))
  
  Au.in[Au.in==" "] <- NA
  Aus <- unique(na.omit(Au.in$AU_STHD))

  v12DMn = NULL
  v13DMn = NULL
  v16DMn = NULL
  v18DMn = NULL
  v20DMn = NULL 
  v22DMn = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL

  for (i in 1:length(Aus)) 
    { 
      SumAU <- SumSummAU[SumSummAU$AU_STHD == Aus[i] & SumSummAU$STHD_IP == "Steelhead",]
      v12Mn = mean(unlist(SumAU$PctDays12))
      v13Mn = mean(unlist(SumAU$PctDays13))
      v16Mn = mean(unlist(SumAU$PctDays16))
      v18Mn = mean(unlist(SumAU$PctDays18))
      v20Mn = mean(unlist(SumAU$PctDays20))
      v22Mn = mean(unlist(SumAU$PctDays22))
      MnMax = mean(unlist(SumAU$MnMx))
      SeMax = MnMax/sqrt(length(SumAU$MnMx))
      
      v12DMn = append(v12DMn, v12Mn)
      v13DMn = append(v13DMn, v13Mn)
      v16DMn = append(v16DMn, v16Mn)
      v18DMn = append(v18DMn, v18Mn)
      v20DMn = append(v20DMn, v20Mn)
      v22DMn = append(v22DMn, v22Mn)
      MnMxDMn = append(MnMxDMn, MnMax)
      SeMxDMn = append(SeMxDMn, SeMax)
      AU_Code = append(AU_Code, Aus[i])
    } 

  SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
  SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)

colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
write.csv(SumAU.out, file = paste0(basin, "_AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, ".csv"), row.names = F)


######
#________Chinook____________________  
######

  popn <- "CHNK"
  varName <- paste0(basin, "_", popn)
  CH_rcas <- unique(na.omit(Au.in[Au.in$CHNK_IP == "Chinook", "RCAID"]))

  Au.in[Au.in==" "] <- NA
  Aus <- unique(na.omit(Au.in$AU_CHIN))
  Aus <- Aus[-1]
  

  v12DMn = NULL
  v13DMn = NULL
  v16DMn = NULL
  v18DMn = NULL
  v20DMn = NULL 
  v22DMn = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL
  
    for (i in 1:length(Aus)) 
      { 
        SumAU <- SumSummAU[SumSummAU$AU_CHIN == Aus[i] & SumSummAU$CHNK_IP == "Chinook",]
        v12Mn = mean(unlist(SumAU$PctDays12))
        v13Mn = mean(unlist(SumAU$PctDays13))
        v16Mn = mean(unlist(SumAU$PctDays16))
        v18Mn = mean(unlist(SumAU$PctDays18))
        v20Mn = mean(unlist(SumAU$PctDays20))
        v22Mn = mean(unlist(SumAU$PctDays22))
        MnMax = mean(unlist(SumAU$MnMx))
        SeMax = MnMax/sqrt(length(SumAU$MnMx))
        
        v12DMn = append(v12DMn, v12Mn)
        v13DMn = append(v13DMn, v13Mn)
        v16DMn = append(v16DMn, v16Mn)
        v18DMn = append(v18DMn, v18Mn)
        v20DMn = append(v20DMn, v20Mn)
        v22DMn = append(v22DMn, v22Mn)
        MnMxDMn = append(MnMxDMn, MnMax)
        SeMxDMn = append(SeMxDMn, SeMax)
        AU_Code = append(AU_Code, Aus[i])
      } 
  
  SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
  SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)
  
  colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
  write.csv(SumAU.out, file = paste0(basin, "_AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, ".csv"), row.names = F)


#######
#_________________pie charts_____________
#######
`

library(plotrix)

basin <- "Aso"
longBasin <- "Asotin"
yrPath <- "13"
yearPath <- "2013"
subDir <- "EP_temp/"
dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
setwd(paste0(mainPath, subDir, longBasin, "/", yearPath))
graphPath <- "D:/OneDrive/work/research/CHaMP/graphics/EP/"

popn <- "CHNK"


  SumAU.out <- read.csv(paste0(basin, "_AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, ".csv"), stringsAsFactors = FALSE)
  
  Aus <- SumAU.out$AU_Code
  metrics <- c("18", "20", "22")

  for (l in metrics)
    {
    
    varName <- paste0("v", l, "DMax")

      for (i in 1:length(Aus)) 
        { 
          ExVal <- SumAU.out[SumAU.out$AU_Code == Aus[i], varName]*100
          value <- c(ExVal, 100-ExVal)
          name <- as.character(Aus[i])
          
          filename <- paste0(graphPath, longBasin, "/", yearPath, "/PctD", l, "/", name, "_", l, popn, ".png", sep="")
          png(filename=filename)
          
          
          if(ExVal >= 10)
          {
            cols <- c("red2", "gainsboro")
            pie(value, col=cols, cex.main=3.0, main=name)
            bisect.angles <- floating.pie(0,0,value, col=cols,)
            pie.labels(0,0,bisect.angles,radius=0.4, c(paste0(value[1],"%")), cex=3, font=2, main=name)
            
          } else if (ExVal < 10 & ExVal > 0){
            cols <- c("red2", "gainsboro")
            pie(value, col=cols, cex.main=3.0, main=name)
            bisect.angles <- floating.pie(0,0,value, col=cols,)
            pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=3, font=2, main=name)
          } else {
            cols <- c("red2", "gainsboro")
            pie(value, clockwise = TRUE, col=cols, cex.main=3.0, main=name)
            bisect.angles <- floating.pie(0,0,value, col=c("gainsboro"), )
            pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=3, font=2, main=name)
          }
          
          dev.off()
        }
    }



