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
  out_Preds <- read.dbf("LST12_USal_interp.dbf") #use the appropriate read statemUSal

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


write.dbf(rca_zonal, file = "LST12_YF_RCA.dbf")





########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/")
  
  LST.in <- read.dbf("LST12_YF_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5))
  colnames(LST.in)[1:46] <- newnamesnum;
  
  elev.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/YF_rca_elev.csv")
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCA_ID")
  
  coeffs.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/2012/All_data_2012_mod_coeffs_Mx.csv")
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
  
  names.out <- sprintf("Tmx_12_%03d", namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]
  
  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/2012/")
  write.dbf(LogPred.out, file = "predt2012_YF_8D_Max_summer.dbf") 

#________________________________________________________
# summary for the summer max EPs
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/"
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

  Au.in <- read.csv("YF_RCA_AUs_IP.csv")
  ST_rcas <- Au.in[Au.in$STHD_IP == "Steelhead",]
  CH_rcas <- Au.in[Au.in$CHNK_IP == "Chinook",]
  
  setwd(paste0(mainPath,yearPath))
  
  Max.in <- read.dbf("predt2012_YF_8D_Max_summer.dbf") 
  
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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2012","Pct13_2012", "Pct16_2012", "Pct18_2012", "Pct20_2012", "Pct22_2012",  "MxMx_2012", "sdMn_2012", "MnMx_2012")


write.dbf(SumSumm.out, file = "YF_2012_21Jul_31Aug_max_summary_ALL.dbf")
write.csv(SumSumm.out, file = "YF_2012_21Jul_31Aug_max_summary_ALL.csv", row.names = F)


# summary by AU -----------------------------------------------------------


  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "RCA_ID")
  colnames(SumSummAU) <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx", "AU_STHD", "AU_CHNK", "STHD_IP", "CHNK_IP")

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
write.csv(SumAU.out, file = "YF_AU_2012_21Jul_31Aug_ExPcnt_summary_ALL.csv", row.names = F)


#______________________
# went back and calc-ed the model metrics for sping & fall together:

# Used the Usal 2012 metrics for this, since there was no new data