############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for a summer max temps MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# 
# Edited Aug 2014 to add the PRESS stastic output
# Should look at PRESS output to screen to chose model - only that model will output resids and preds.

#Edited Oct 2015 for EP watershed summer output

          
##############################################################################################
# This section reads in the 1km LST data for a year, uses a 4th order polynomial to fill in Julian day 1 & 365 (if they are missing)
# then fills any remaining gaps across the year at each pixel with a linear interpolation.
##############################################################################################


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)
library(lubridate)
library(zoo)
library(car)
library(gvlma)

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/shapes/")

  LST.in <- read.dbf("UGR_12_LST.dbf")
  PointID <- LST.in[,1]
  LST.in <- LST.in[,2:47]
  LST.in[LST.in<0.1] = NA
  
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
  
  plot(1:46, LST.out[1, 1:46])
  
  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/")
  write.dbf(LST.out, file = "LST12_UGR_interp.dbf")

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/")
  out_Preds <- read.dbf("LST12_UGR_interp.dbf") #use the appropriate read statemUGR
  
  weights <- read.csv("UGR_rca_area_wgts.csv")
  weights[weights$area_wgts<0.01, "area_wgts"] = 0.01
  
  rcas <- unique(unlist(weights$RCAID))
  rca_zonal <- matrix(ncol=46, nrow=length(rcas))


  l <- 1
  for(i in rcas)	#assumes the RCAs are consecutively numbered starting with "1"
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

plot(1:46, rca_zonal[1, 1:46])

write.dbf(rca_zonal, file = "LST12_UGR_RCA.dbf")



#################################################################
#Logger prediction modeling part
#Generates models for both LST-only and LST + Julian Day datasets
#################################################################
  
  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/")

  ID.in <- read.csv("UGR_sites_elev.csv")
     
  LST.in <- read.dbf("LST12_UGR_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    
  colnames(LST.in)[1:46] <- newnamesnum;

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/"
  yearPath <- "2012"

  setwd(paste0(mainPath,yearPath))
  Log.in <- read.csv("UGR_2012_8DMax_logger.csv")
  Log.in$SiteName <- as.character(Log.in$SiteName)
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)
  
 
  LST.Log.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteName  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$RCAID[ID.in$SiteName == i]
      Elev <- ID.in$Elev_M[ID.in$SiteName == i]
      
      LST.site <- matrix(ncol=3, nrow=46)
      LST.site[,1] <- newnamesnum
      LST.site[,2] <- unlist(LST.in[RCAID,1:46])
      LST.site <- data.frame(LST.site)
      colnames(LST.site) <- c("JulDay", "LST", "Elev")
      LST.site[3] <- Elev
      LST.Log.site <- merge(LST.site, Log.site, by.x = "JulDay", by.y = "JulDay", all.x=TRUE, all.y = FALSE)
      
      LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
    }
    

    ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
    NoNA.xyz <- LST.Log.out[ind,]

    NoNA.xyz <- NoNA.xyz[,c(8, 2, 1, 3, 5)]
    colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")

    write.csv(x=NoNA.xyz, file="UGR_2012_8Day_model_data.csv", row.names = FALSE)
    NoNA.xyz$y <- as.numeric(NoNA.xyz$y)
    NoNA.xyz <- orderBy(~z, NoNA.xyz)

#################################
# to suss out summer stuff
###############################

  summer <- subset(NoNA.xyz, z > 181 & z < 258)
  SiteID <- as.matrix(unique(summer$SiteName))
  y <- summer$y
  x <- summer$x
  z <- summer$z
  e <- summer$e
  plot (x,y)

  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y<0] = 0.0
  
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
  RSS <- sum(pressstat_sum$residuals^2)
  TSS <- sum((y - mean(y))^2)
  MSS <- TSS - RSS
  p2 <- MSS/TSS

  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "summer"
  pred.out[,5] <- "2012"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_UGR_2012_summer_only.csv",sep = ",", col.names=F)  

  plot(pred.out[,1], pred.out[,2])
  abline(0,1)

################################
# full year
###################################

  y <- NoNA.xyz$y
  x <- NoNA.xyz$x
  z <- NoNA.xyz$z
  e <- NoNA.xyz$e
  plot(x, y)

  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y<0] = 0.0
  
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
  RSS <- sum(pressstat_sum$residuals^2)
  TSS <- sum((y - mean(y))^2)
  MSS <- TSS - RSS
  p2 <- MSS/TSS
  
  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "full year"
  pred.out[,5] <- "2012"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_UGR_2012_full_year.csv",sep = ",", col.names=F)  

  summer_pred <- subset(pred.out, z > 181 & z < 258)
  points(summer_pred[, 1], summer_pred[,2], pch = 16, col = "green")

#####################################
# spring/fall
######################################

    coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bJul=numeric(2), bElev=numeric(2))
    metrics_out <- data.frame(PRESS=numeric(2), p2=numeric(2), r2=numeric(2), RMSEP=numeric(2), RMSE=numeric(2))
    rownames(metrics_out) <- c("Spring", "Fall")
    rownames(coeffs_out) <- c("Spring", "Fall")

    maxrow <- which.max(NoNA.xyz$y)
    data.sp <- NoNA.xyz[1:maxrow,]
    data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]
    
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
    pred.y[pred.y<0] = 0.0
    
    plot(pred.y, y)
    gvmodel <- gvlma(mod)
    summary(gvmodel)
    qqPlot(mod, main="QQ Plot")
    spreadLevelPlot(mod)
    plot(mod$residuals, x)
        
    pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
    RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
    RSS <- sum(pressstat_sum$residuals^2)
    TSS <- sum((y - mean(y))^2)
    MSS <- TSS - RSS
    p2 <- MSS/TSS

    pred.out <- matrix(nrow=length(pred.y), ncol=5)
    pred.out[,1] <- y
    pred.out[,2] <- pred.y
    pred.out[,3] <- z
    pred.out[,4] <- "spring"
    pred.out[,5] <- "2012"
    colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
    write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_UGR_2012_sp_fall.csv",sep = ",", col.names=F)  

    sp_pred <- subset(pred.out, z > 181 & z < 258)
    points(sp_pred[, 1], sp_pred[,2], pch = 16, col = "red")
    plot(sp_pred[, 1], sp_pred[,2])
    abline(0,1)
    
    coeffs_out[1,1] <- coeffs[1,1]
    coeffs_out[1,2] <- coeffs[2,1]
    coeffs_out[1,3] <- coeffs[3,1]
    coeffs_out[1,4] <- coeffs[4,1]
    coeffs_out[1,5] <- coeffs[5,1]
    
    metrics_out[1,1] <- pressstat_sum$stat
    metrics_out[1,2] <- pressstat_sum$P.square
    metrics_out[1,3] <- sum_mod$adj.r.squared
    metrics_out[1,4] <- RMSEP
    metrics_out[1,5] <- sum_mod$sigma
    

  
    y <- data.fall$y
    x <- data.fall$x
    z <- data.fall$z
    e <- data.fall$e
    mod <- lm(y ~ x + I(x^2) + z + e)
    sum_mod <- summary(mod)
    coeffs <- as.matrix(coefficients(mod))
    pred.y <- predict(mod)
    pred.y[pred.y<0] = 0.0

    plot(pred.y, y)
    gvmodel <- gvlma(mod)
    summary(gvmodel)
    qqPlot(mod, main="QQ Plot")
    spreadLevelPlot(mod)
    plot(mod$residuals, x)
    
    
    pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
    RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
    RSS <- sum(pressstat_sum$residuals^2)
    TSS <- sum((y - mean(y))^2)
    MSS <- TSS - RSS
    p2 <- MSS/TSS

    pred.out <- matrix(nrow=length(pred.y), ncol=5)
    pred.out[,1] <- y
    pred.out[,2] <- pred.y
    pred.out[,3] <- z
    pred.out[,4] <- "fall"
    pred.out[,5] <- "2012"
    
    write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_UGR_2012_sp_fall.csv",sep = ",", col.names=F)  

    fall_pred <- subset(pred.out, z > 181 & z < 258)
    origY <- fall_pred[,1]
    predY <- fall_pred[,2]
    points(fall_pred[, 1], fall_pred[,2], pch = 16, col = "purple")
    
           
    coeffs_out[2,1] <- coeffs[1,1]
    coeffs_out[2,2] <- coeffs[2,1]
    coeffs_out[2,3] <- coeffs[3,1]
    coeffs_out[2,4] <- coeffs[4,1]
    coeffs_out[2,5] <- coeffs[5,1]
    
    metrics_out[2,1] <- pressstat_sum$stat
    metrics_out[2,2] <- p2
    metrics_out[2,3] <- sum_mod$adj.r.squared
    metrics_out[2,4] <- RMSEP
    metrics_out[2,5] <- sum_mod$sigma

write.table(x=coeffs_out, append=F,row.names=T, file = "Summer_data_2012_mod_coeffs_Mx.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = "Summer_data_2012_mod_metrics_Mx.csv", sep = ",", col.names=T)



########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################


  
  LST.in <- read.dbf("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/LST12_UGR_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5))
  colnames(LST.in)[1:46] <- newnamesnum;
  
  elev.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/UGR_rca_elev.csv")
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCAID")
  
  coeffs.in <- read.csv("Summer_data_2012_mod_coeffs_Mx.csv")
  LogPred.out <- LST.elev[,c(25:36,48,1)]
  LogPred.out[,1:12] <- 0
  rcas <- unique(elev.in$RCAID)
  LST.sum <- LST.elev[,c(25:36, 48,1)]

#     for (i in 1:length(rcas))  
#      {
#         x <- unlist(LST.sum[i,])
#         maxrow <- as.numeric(which.max(x[1:12])) #either specify or let be dynamic
#         day <- as.numeric(colnames(LST.sum)[maxrow])
#         
#         j <- 185
#         for (l in 1:maxrow)
#         {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
#          j <- j + 8}
#         k <- day
#         for (l in maxrow:12)     
#         {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[13] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
#          k <- k + 8}
#         LogPred.out[i,1:12] <- x [1:12] 
#       }

-
      for (i in 1:length(rcas))  
        {
          x <- unlist(LST.sum[i,])
          maxrow <- as.numeric(which.max(x[1:12])) #either specify or let be dynamic
          day <- as.numeric(colnames(LST.sum)[maxrow])
          
          j <- 185
          for (l in 1:12)
          {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
           j <- j + 8}
          LogPred.out[i,1:12] <- x [1:12] 
        }  

    LogPred.out <- as.data.frame(LogPred.out)
    
  LogPred.out[LogPred.out< -0.5] = -0.5
  LogPred.out$Basin_RCA <- paste0("UGR_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:12]))
  
  names.out <- sprintf("Tmx_12_%03d", namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]
  
  plot(1:12, LogPred.out[1,1:12])
  points(1:12, LogPred.out[100,1:12], pch = 16, col = "blue")

  write.dbf(LogPred.out, file = "predt2012_UGR_8D_Max_summer.dbf") 

#________________________________________________________
# summary for the summer max EPs
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/"
yearPath <- "2012"


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

  setwd(mainPath)
  
  Au.in <- read.csv("UGR_RCA_AUs_IP.csv", stringsAsFactors=F)
 
  
  setwd(paste0(mainPath,yearPath))
  
  Max.in <- read.dbf("predt2012_UGR_8D_Max_summer.dbf") 
  
  
  Max.au <- merge(Max.in, Au.in, by.x = "RCAID", by.y = "RCAID")
  
  Max.daily <- Max.au[, c(1, rep(4:9, each = 8), 16:19)]
  colnames(Max.daily)[2:49] <- 201:248


#________________All____________________________________________
  rcas <- unique(Au.in$RCAID)
#__________________________________________________________
  
  SumSumm.out <- data.frame ("RCAID" = rcas)


# count of days in exceedence -----------------------------------------


    for (i in 1:length(rcas)) 
      { 
        l <- rcas[i]
        MaxRCA <- Max.daily[Max.daily$RCAID == l,] #grab days for one RCA 
        DaysAbove12 <- length(which(MaxRCA[2:43]> 12)) #finds how many day in the 20July-31Aug window exceed threshold 
        DaysAbove13 <- length(which(MaxRCA[2:43]> 13))
        DaysAbove16 <- length(which(MaxRCA[2:43]> 16))
        DaysAbove18 <- length(which(MaxRCA[2:43]> 18))
        DaysAbove20 <- length(which(MaxRCA[2:43]> 20))
        DaysAbove22 <- length(which(MaxRCA[2:43]> 22))
        MaxMax <- max(MaxRCA[2:43])
        SDMax <- sd(MaxRCA[2:43])
        MeanMax <- mean(unlist(MaxRCA[2:43]))
        SumSumm.out$PctDays12[i] <- DaysAbove12/42
        SumSumm.out$PctDays13[i] <- DaysAbove13/42
        SumSumm.out$PctDays16[i] <- DaysAbove16/42
        SumSumm.out$PctDays18[i] <- DaysAbove18/42
        SumSumm.out$PctDays20[i] <- DaysAbove20/42
        SumSumm.out$PctDays22[i] <- DaysAbove22/42
        SumSumm.out$MxMx[i] <- MaxMax 
        SumSumm.out$SdMn[i] <- SDMax
        SumSumm.out$MnMx[i] <- MeanMax
      } 


SumSumm.out[,2:10] <- round(SumSumm.out[,2:10], digits = 2)

colnames(SumSumm.out) <- c("RCAID", "Pct12_2012","Pct13_2012", "Pct16_2012", "Pct18_2012", "Pct20_2012", "Pct22_2012",  "MxMx_2012", "sdMn_2012", "MnMx_2012")


write.dbf(SumSumm.out, file = "UGR_2012_21Jul_31Aug_max_summary_All.dbf")
#write.csv(SumSumm.out, file = "UGR_2012_21Jul_31Aug_max_summary_All.csv", row.names = F)


# summary by AU -----------------------------------------------------------

SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "RCAID")
colnames(SumSummAU)[1:10] <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx")

#_______________Chinook__________________________________________  
# summmary by AU 
#______________________________________

  CH_Aus <- as.character(unique(Au.in$AU_CHIN))    

  v12DMn = NULL
  v13DMn = NULL
  v16DMn = NULL
  v18DMn = NULL
  v20DMn = NULL 
  v22DMn = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL


  for (i in 1:length(CH_Aus)) 
    { 
      SumAU <- SumSummAU[SumSummAU$AU_CHIN == CH_Aus[i],]
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
      AU_Code = append(AU_Code, CH_Aus[i])
    } 

SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)

colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
write.csv(SumAU.out, file = "UGR_AU_2012_21Jul_31Aug_ExPcnt_summary_CHNK_All.csv", row.names = F)

#________________________________
# summmary of IP reaches only by AU for Chinook
#______________________________________


  v12DMn = NULL
  v13DMn = NULL
  v16DMn = NULL
  v18DMn = NULL
  v20DMn = NULL 
  v22DMn = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL


    for (i in 1:length(CH_Aus)) 
      { 
        SumAU <- SumSummAU[SumSummAU$AU_CHIN == CH_Aus[i] & SumSummAU$CHNK_IP == "Chinook",]
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
        AU_Code = append(AU_Code, CH_Aus[i])
      } 

SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)

colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
write.csv(SumAU.out, file = "UGR_AU_2012_21Jul_31Aug_ExPcnt_summary_CHNK_IP.csv", row.names = F)

#_______________Steelhead___________________________________________ 
# summmary by AU for Steelhead
#______________________________________

ST_Aus <- as.character(unique(ST_Au_rcas$AU_STHD))

  v12DMn = NULL
  v13DMn = NULL
  v16DMn = NULL
  v18DMn = NULL
  v20DMn = NULL 
  v22DMn = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL
  
  
    for (i in 1:length(ST_Aus)) 
      { 
        SumAU <- SumSummAU[SumSummAU$AU_STHD == ST_Aus[i],]
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
        AU_Code = append(AU_Code, ST_Aus[i])
      } 

SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)

colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
write.csv(SumAU.out, file = "UGR_AU_2012_21Jul_31Aug_ExPcnt_summary_STHD_All.csv", row.names = F)


#________________________________
# summmary of IP reaches only by AU for Steelhead
#______________________________________

  v12DMn = NULL
  v13DMn = NULL
  v16DMn = NULL
  v18DMn = NULL
  v20DMn = NULL 
  v22DMn = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL

    for (i in 1:length(ST_Aus)) 
      { 
        SumAU <- SumSummAU[SumSummAU$AU_STHD == ST_Aus[i] & SumSummAU$STHD_IP == "Steelhead",]
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
        AU_Code = append(AU_Code, ST_Aus[i])
      } 

SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)

colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
write.csv(SumAU.out, file = "UGR_AU_2012_21Jul_31Aug_ExPcnt_summary_STHD_IP.csv", row.names = F)


#________________________________________________________
# summary for the summer max EPs linear interpolation
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/UGR/"
yearPath <- "2012"


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

setwd(mainPath)

  Au.in <- read.csv("UGR_RCA_AUs_IP.csv")
  ST_rcas <- read.csv("UGR_RCA_AUs_STHD.csv")
  CH_rcas <- read.csv("UGR_RCA_AUs_CHNK.csv")

  colnames(ST_rcas)[2] <- "AU_STHD"
  
  setwd(paste0(mainPath,yearPath))

  Max.in <- read.dbf("predt2012_UGR_8D_Max_summer.dbf") 

  Max.in$C <- as.character(NA)
  Max.rep <- Max.in[,c(3,16,4,16,5,16,6,16,7,16,8,16,9)]
  toFill <- t(Max.rep)
  filled <- t(interpNA(toFill, method="linear"))
  dim(filled)
  
  toFill <- cbind(filled, Max.in[,16])
  toFill.rep <- toFill[,c(1,14,2,14,3,14,4,14,5,14,6,14,7,14,8,14,9,14,10,14,11,14,12,14,13)]
  ttoFill <- t(toFill.rep)
  filled <- t(interpNA(ttoFill, method="linear"))
  dim(filled)
  
  toFill <- cbind(filled, Max.in[,16])
  toFill.rep <- toFill[,c(1,26,2,26,3,26,4,26,5,26,6,26,7,26,8,26,9,26,10,26,11,26,12,26,13,26,14,26,15,26,16,26,17,26,18,26,19,26,20,26,21,26,22,26,23,26,24,26,25)]
  ttoFill <- t(toFill.rep)
  filled <- t(interpNA(ttoFill, method="linear"))
  dim(filled)
  
  class(filled) <- "numeric"
  filled.out <- as.data.frame(filled, stringsAsFactors=FALSE)
  
  filled.out$RCAID <- Max.in$RCAID
  colnames(filled.out)[1:49] <- 201:249

  

  Max.daily <- filled.out

  #____________________________  
  Aus <- as.character(unique(ST_rcas$AU_STHD))
  rcas <- unique(ST_rcas$RCAID)
  #_______________________________
  Aus <- as.character(unique(CH_rcas$AU_CHNK))
  rcas <- unique(CH_rcas$RCAID)
  #________________________________  
  rcas <- unique(Au.in$RCAID)

SumSumm.out <- data.frame ("RCAID" = rcas)


# count of days in exceedence -----------------------------------------


  for (i in 1:length(rcas)) 
    { 
      l <- rcas[i]
      MaxRCA <- Max.daily[Max.daily$RCAID == l,] #grab days for one RCA 
      DaysAbove12 <- length(which(MaxRCA[1:43]> 12)) #finds how many day in the 20July-31Aug window exceed threshold 
      DaysAbove13 <- length(which(MaxRCA[1:43]> 13))
      DaysAbove16 <- length(which(MaxRCA[1:43]> 16))
      DaysAbove18 <- length(which(MaxRCA[1:43]> 18))
      DaysAbove20 <- length(which(MaxRCA[1:43]> 20))
      DaysAbove22 <- length(which(MaxRCA[1:43]> 22))
      MaxMax <- max(MaxRCA[1:43])
      SDMax <- sd(MaxRCA[1:43])
      MeanMax <- mean(unlist(MaxRCA[1:43]))
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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2012","Pct13_2012", "Pct16_2012", "Pct18_2012", "Pct20_2012", "Pct22_2012",  "MxMx_2012", "sdMn_2012", "MnMx_2012")

write.dbf(SumSumm.out, file = "UGR_2012_21Jul_31Aug_max_summary_STHD_redo.dbf")
write.csv(SumSumm.out, file = "UGR_2012_21Jul_31Aug_max_summary_STHD_redo.csv", row.names = F)


# summary by AU -----------------------------------------------------------


SumSummAU <- merge(SumSumm.out, CH_rcas, by.x = "RCAID", by.y = "RCAID", all.x=TRUE, all.y=FALSE)
colnames(SumSummAU)[11] <- "AU_Code"
SumSummAU$AU_Code <- as.character(SumSummAU$AU_Code)

colnames(SumSummAU) <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx", "AU_Code")

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
      SumAU <- SumSummAU[SumSummAU$AU_Code == Aus[i],]
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
      AU_Code = append(AU_Code, as.character(Aus[i]))
    } 

SumAU.out = data.frame(v12DMn, v13DMn, v16DMn, v18DMn, v20DMn, v22DMn, MnMxDMn, SeMxDMn, AU_Code)
SumAU.out[,1:8] <- round(SumAU.out[,1:8], digits = 2)

colnames(SumAU.out) <- c("v12DMax", "v13DMax", "v16DMax", "v18DMax", "v20DMax", "v22DMax", "MnMax", "SeMnMx", "AU_Code")
write.csv(SumAU.out, file = "UGR_AU_2012_21Jul_31Aug_ExPcnt_summary_CHNK_redo.csv", row.names = F)


#_________________pie charts_____________


library(plotrix)

  SumAU.out$AU_Code <- as.character(SumAU.out$AU_Code)
  

  for (i in 1:length(Aus)) 
    { 
      ExVal <- SumAU.out[SumAU.out$AU_Code == Aus[i], 'v18DMax']*100
      value <- c(ExVal, 100-ExVal)
      name <- as.character(Aus[i])
      
      filename <- paste0("D:/OneDrive/work/research/CHaMP/graphics/EP/UGR/2012/PctD18/CHNK/", name, "_18CHNK.png", sep="")
      png(filename=filename)
      
      
      if(ExVal >= 10)
        {
          cols <- c("red2", "gainsboro")
          pie(value, col=cols, cex.main=2.0, main=name)
          bisect.angles <- floating.pie(0,0,value, col=cols)
          pie.labels(0,0,bisect.angles,radius=0.4, c(paste0(value[1],"%")), cex=2.5, font=2, main=name)
          
        } else if (ExVal < 10 & ExVal > 0){
          cols <- c("red2", "gainsboro")
          pie(value, col=cols, cex.main=2.0, main=name)
          bisect.angles <- floating.pie(0,0,value, col=cols)
          pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=2.5, font=2, main=name)
        } else {
          cols <- c("red2", "gainsboro")
          pie(value, clockwise = TRUE, col=cols, cex.main=2.0, main=name)
          bisect.angles <- floating.pie(0,0,value, col=c("gainsboro"), )
          pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=2.5, font=2, main=name)
        }
      
      dev.off()
    }
