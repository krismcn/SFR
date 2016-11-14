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


####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/shapes/")
  out_Preds <- read.dbf("LST13_USal_interp.dbf") #use the appropriate read statemUSal

  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/")
  weights <- read.csv("YF_area_wgts.csv")
  weights[weights$area_wgts<0.01, "area_wgt"] = 0.01
  
  rcas <- unique(unlist(weights$RCA_ID))
  rca_zonal <- matrix(ncol=46, nrow=length(rcas))


  l <- 1
  for(i in rcas)	#assumes the RCAs are consecutively numbered starting with "1"
      {
      pixels <- weights[weights$RCA_ID == i, "GRIDCODE"]
      wgts <- weights[weights$RCA_ID == i, "area_wgt"]
      
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


write.dbf(rca_zonal, file = "LST13_YF_RCA.dbf")

#################################################################
#Logger prediction modeling part
#Generates models for both LST-only and LST + Julian Day datasets
#################################################################

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/")

  ID.in <- read.csv("YF_sites_elev.csv")
  
  LST.in <- read.dbf("LST13_YF_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    
  colnames(LST.in)[1:46] <- newnamesnum;
  
  
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/"
  yearPath <- "2013"

  setwd(paste0(mainPath,yearPath))
  Log.in <- read.csv("YF_2013_logger.csv")
  Log.in$SiteName <- as.character(Log.in$SiteName)
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)


  LST.Log.out <- data.frame (mup = NULL)

    for (i in SiteID) 
      { 
        Log.site <- Log.in[Log.in$SiteName  == i,]
        Log.site <- as.data.frame(Log.site)
        
        
        RCAID <- ID.in$RCA_ID[ID.in$SiteName == i]
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
  plot(NoNA.xyz$x, NoNA.xyz$y)

  write.csv(x=NoNA.xyz, file="YF_2013_8Day_model_data.csv", row.names = FALSE)
  NoNA.xyz$y <- as.numeric(NoNA.xyz$y)
  NoNA.xyz <- orderBy(~z, NoNA.xyz)

####spring#######

  maxrow <- which.max(NoNA.xyz$y)
  y <- NoNA.xyz$y[1:maxrow]
  x <- NoNA.xyz$x[1:maxrow]
  z <- NoNA.xyz$z[1:maxrow]
  e <- NoNA.xyz$e[1:maxrow]
  plot(x,y)

####fall#######

  maxrow <- which.max(NoNA.xyz$y)
  y <- NoNA.xyz$y[maxrow:nrow(NoNA.xyz)]
  x <- NoNA.xyz$x[maxrow:nrow(NoNA.xyz)]
  z <- NoNA.xyz$z[maxrow:nrow(NoNA.xyz)]
  e <- NoNA.xyz$e[maxrow:nrow(NoNA.xyz)]
  plot(x,y)

########################
  summer <- subset(NoNA.xyz, z > 181 & z < 258)
  SiteID <- as.matrix(unique(summer$SiteName))
  SiteID <- as.matrix(as.character(SiteID))
 
  
  y <- summer$y
  x <- summer$x
  z <- summer$z
  e <- summer$e
  plot(x, y)

####All#######################    
  y <- NoNA.xyz$y
  x <- NoNA.xyz$x
  z <- NoNA.xyz$z
  e <- NoNA.xyz$e
  plot(x, y)

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
    mod <- lm(y ~ x + I(x^2) + z + e)
    sum_mod <- summary(mod)
    coeffs <- as.matrix(coefficients(mod))
    pred.y <- predict(mod)
    pred.y[pred.y<0] = 0.0
    
    
    
    pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
    RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
    
    resids <- pressstat_sum$residuals
    
    pred.out <- matrix(nrow=length(pred.y), ncol=5)
    pred.out[,1] <- y
    pred.out[,2] <- pred.y
    pred.out[,3] <- z
    pred.out[,4] <- "all"
    pred.out[,5] <- "2013"
    colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
    write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_YF_2013_yr.csv",sep = ",", col.names=F)  
    
    
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
  
    all_resids <- resids
    
    RSS <- sum(all_resids^2)
    TSS <- sum((y - mean(y))^2)
    MSS <- TSS - RSS
    p2 <- MSS/TSS
    RMSEP <- sqrt(mean(all_resids^2))

    metrics_out[1,2] <- p2
    coeffs_out[2,] <- coeffs_out[1,]

write.table(x=coeffs_out, append=F,row.names=T, file = "All_data_2013_mod_coeffs_Mx.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = "All_data_2011_mod_metrics_Mx.csv", sep = ",", col.names=T)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/")
  
  LST.in <- read.dbf("LST13_YF_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5))
  colnames(LST.in)[1:46] <- newnamesnum;
  
  elev.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/YF_rca_elev.csv")
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCA_ID")
  
  coeffs.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/2013/All_data_2013_mod_coeffs_Mx.csv")
  LogPred.out <- LST.elev[,c(25:36,48,1)]
  LogPred.out[,1:12] <- 0
  rcas <- unique(elev.in$RCA_ID)
  LST.sum <- LST.elev[,c(25:36, 48,1)]

    for (i in 1:length(rcas))  
     {
        x <- unlist(LST.sum[i,])
        maxrow <- as.numeric(which.max(x[1:12])) #either specify or let be dynamic
        day <- as.numeric(colnames(LST.sum)[maxrow])
        
        j <- 185
        for (l in 1:maxrow)
        {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
         j <- j + 8}
        k <- day
        for (l in maxrow:12)     
        {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[13] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
         k <- k + 8}
        LogPred.out[i,1:12] <- x [1:12] 
      }
    
    LogPred.out <- as.data.frame(LogPred.out)
    
  LogPred.out[LogPred.out< -0.5] = -0.5
  LogPred.out$Basin_RCA <- paste0("YF_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:12]))
  
  names.out <- sprintf("Tmx_13_%03d", namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]
  
  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/2013/")
  write.dbf(LogPred.out, file = "predt2013_YF_8D_Max_summer.dbf") 

#________________________________________________________
# summary for the summer max EPs
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/"
yearPath <- "2013"


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

  setwd(mainPath)
  
  Au.in <- read.csv("YF_RCA_AUs_IP.csv")
  ST_rcas <- Au.in[Au.in$STHD_IP == "Steelhead",]
  CH_rcas <- Au.in[Au.in$CHNK_IP == "Chinook",]
  
  setwd(paste0(mainPath,yearPath))
  
  Max.in <- read.dbf("predt2013_YF_8D_Max_summer.dbf") 
  
  Max.au <- merge(Max.in, Au.in, by.x = "RCAID", by.y = "RCA_ID")
  
  Max.au$AU_STHD <- as.character(Max.au$AU_STHD)
  Max.au$AU_CHNK <- as.character(Max.au$AU_CHIN)
  

  Max.daily <- Max.in[, c(rep(3:8, each = 8), 13:14)]
  colnames(Max.daily)[1:48] <- 201:248
#____________________________  
  Aus <- unique(Max.au$AU_STHD)
  rcas <- unique(ST_rcas$RCA_ID)
#_______________________________
  Aus <- unique(Max.au$AU_CHNK)
  rcas <- unique(CH_rcas$RCA_ID)
#________________________________  
  rcas <- unique(Max.au$RCAID)

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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2013","Pct13_2013", "Pct16_2013", "Pct18_2013", "Pct20_2013", "Pct22_2013",  "MxMx_2013", "sdMn_2013", "MnMx_2013")

write.dbf(SumSumm.out, file = "YF_2013_21Jul_31Aug_max_summary_ALL.dbf")
write.csv(SumSumm.out, file = "YF_2013_21Jul_31Aug_max_summary_ALL.csv", row.names = F)



# summary by AU -----------------------------------------------------------


  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "RCA_ID")
  SumSummAU <- SumSummAU[,1:12]
  colnames(SumSummAU) <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx", "AU_STHD", "AU_CHNK")
 
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
      SumAU <- SumSummAU[SumSummAU$AU_CHNK == Aus[i],]
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
write.csv(SumAU.out, file = "YF_AU_2013_21Jul_31Aug_ExPcnt_summary_All.csv", row.names = F)

#______________________
# went back and calc-ed the model metrics for sping & fall together:

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/2013/")

  Log.in <- read.csv("YF_2013_8Day_model_data.csv")
  NoNA.xyz <- Log.in
  NoNA.xyz$SiteName <- as.character(NoNA.xyz$SiteName)
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:maxrow,]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]
  y <- data.sp$y
  x <- data.sp$x
  z <- data.sp$z
  e <- data.sp$e
  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  
  resids <- pressstat_sum$residuals
  y <- data.fall$y
  x <- data.fall$x
  z <- data.fall$z
  e <- data.fall$e
  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  
  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  
  all_resids <- c(resids, pressstat_sum$residuals)
  y <- c(data.sp$y, data.fall$y)
  RSS <- sum(all_resids^2)
  TSS <- sum((y - mean(y))^2)
  MSS <- TSS - RSS
  p2 <- MSS/TSS
  RMSEP <- sqrt(mean(all_resids^2))

#___________________________

#________________________________________________________
# summary for the summer max EPs no reps
# 15Jul-31Aug ------------------------------------------------------------------



  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/"
  yearPath <- "2013"


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

  setwd(mainPath)

  Au.in <- read.csv("YF_RCA_AUs_IP.csv")
  ST_rcas <- Au.in[Au.in$STHD_IP == "Steelhead",]
  CH_rcas <- Au.in[Au.in$CHNK_IP == "Chinook",]
  
  setwd(paste0(mainPath,yearPath))
  
  Max.in <- read.dbf("predt2013_YF_8D_Max_summer.dbf") 
   
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
  
  filled.out$RCAID <- Max.daily$RCAID
  colnames(filled.out)[1:49] <- 201:249

  Max.au <- merge(Max.in, Au.in, by.x = "RCAID", by.y = "RCA_ID")
  
  Max.au$AU_STHD <- as.character(Max.au$AU_STHD)
  Max.au$AU_CHNK <- as.character(Max.au$AU_CHIN)


  Max.daily <- filled.out
  
  #____________________________  
  Aus <- unique(Max.au$AU_STHD)
  rcas <- unique(ST_rcas$RCA_ID)
  #_______________________________
  Aus <- unique(Max.au$AU_CHNK)
  rcas <- unique(CH_rcas$RCA_ID)
  #________________________________  
  rcas <- unique(Max.au$RCAID)
  
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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2013","Pct13_2013", "Pct16_2013", "Pct18_2013", "Pct20_2013", "Pct22_2013",  "MxMx_2013", "sdMn_2013", "MnMx_2013")

write.dbf(SumSumm.out, file = "YF_2013_21Jul_31Aug_max_summary_ALL_redo.dbf")
write.csv(SumSumm.out, file = "YF_2013_21Jul_31Aug_max_summary_ALL_redo.csv", row.names = F)


# summary by AU -----------------------------------------------------------


  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "RCA_ID")
  SumSummAU <- SumSummAU[,1:12]
  colnames(SumSummAU) <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx", "AU_STHD", "AU_CHNK")
  
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
      SumAU <- SumSummAU[SumSummAU$AU_CHNK == Aus[i],]
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
write.csv(SumAU.out, file = "YF_AU_2013_21Jul_31Aug_ExPcnt_summary_CHNK_redo.csv", row.names = F)

