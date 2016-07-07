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



########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/")
  
  LST.in <- read.dbf("LST11_YF_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5))
  colnames(LST.in)[1:46] <- newnamesnum;
  
  elev.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/YF_rca_elev.csv")
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCA_ID")
  
  coeffs.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/2011/All_data_2011_mod_coeffs_Mx.csv")
  LogPred.out <- LST.elev[,c(21:36,48,1)]
  LogPred.out[,1:16] <- 0
  rcas <- unique(elev.in$RCA_ID)
  LST.sum <- LST.elev[,c(21:36, 48,1)]

    for (i in 1:length(rcas))  
     {
        x <- unlist(LST.sum[i,])
        maxrow <- as.numeric(which.max(x[1:16])) #either specify or let be dynamic
        day <- as.numeric(colnames(LST.sum)[maxrow])
        
        j <- 153
        for (l in 1:maxrow)
        {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[17] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
         j <- j + 8}
        k <- day
        for (l in maxrow:16)     
        {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[17] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
         k <- k + 8}
        LogPred.out[i,1:16] <- x [1:16] 
      }
    
    LogPred.out <- as.data.frame(LogPred.out)
    
  LogPred.out[LogPred.out< -0.5] = -0.5
  LogPred.out$Basin_RCA <- paste0("YF_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:16]))
  
  names.out <- sprintf("Tmx_11_%03d", namesnum)
  colnames(LogPred.out)[1:16] <- names.out[1:16]
  
  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/YF/2011/")
  write.dbf(LogPred.out, file = "predt2011_YF_8D_Max_ext_summer_13_coeffs.dbf") 

