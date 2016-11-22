############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for a year using MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent

# Edited Aug 2014 to add the PRESS stastic output
# Should look at PRESS output to screen to chose model - only that model will output resids and preds.


          
##############################################################################################
# This section reads in the 1km LST data 
# then fills any gaps across the year at each pixel with an interpolation.
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

mainPath <- "D:/OneDrive/work/research/CHaMP/"
basin <- "Ent"
longBasin <- "Entiat"
subDir <- "CHaMP_data/"
setwd(paste0(mainPath, subDir, longBasin))

  LST.in <- read.csv("Ent_LST_2013.csv", header = TRUE, stringsAsFactors=FALSE)
  GrPolID <- LST.in[,1]
  LST.in <- LST.in[,2:47]
  LST.in[LST.in<0.1] = NA
  
  LST.names<-colnames(LST.in)
  
  tLST.in <- t(LST.in)
  tLST.out <- na.spline(tLST.in)
  CLST <- tLST.out*0.02-273.15
  tCLST <- t(CLST)
  CLST.out <- as.data.frame(tCLST)
  plot(1:46, CLST.out[10, 1:46])
  CLST.out$GRID_CODE <- GrPolID
  
  colnames(CLST.out)[1:46] <- LST.names
    
  
  write.dbf(CLST.out, file = "LST13_Ent_interp.dbf")

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

out_Preds <- read.dbf("LST13_Ent_interp.dbf") #use the appropriate read statement

  weights <- read.dbf("Ent_rca_area_wgts.dbf")
  weights[weights$area_wgts<0.01, "area_wgts"] = 0.01
  rcas <- unique(unlist(weights$rca_id))
  rca_zonal <- matrix(ncol=46, nrow=length(rcas))
  
    l <- 1
    for(i in rcas)  
      {
        pixels <- weights[weights$rca_id == i, "GRIDCODE"]
        wgts <- weights[weights$rca_id == i, "area_wgt"]
        
        for (j in 1:46)
          {
            daily <- out_Preds[out_Preds$GRID_CODE %in% pixels, j]
            zonal_mn <- weighted.mean(daily, wgts, na.rm = TRUE)
            rca_zonal[l,j] <- zonal_mn
          }
          l <- l+1
      }
  
  colnames(rca_zonal)[1:46] <- colnames(out_Preds)[1:46]
  rca_zonal <- as.data.frame(rca_zonal)
  rca_zonal$RCAID <- rcas
  plot(1:46, rca_zonal[1,1:46])

write.dbf(rca_zonal, file = "LST13_Ent_RCA.dbf")



#################################################################
#Logger prediction modeling part

#################################################################
setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Entiat")

  ID.in <- read.csv("Ent_sites_elev.csv", stringsAsFactors=FALSE)


  LST.in <- read.dbf("LST13_Ent_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- substring (colnames(LST.in)[1:46],3,5)
  newnamesnum <- as.numeric(newnamesnum)
  colnames(LST.in)[1:46] <- newnamesnum
  

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
yearPath <- "2013"
setwd(paste0(mainPath))

  Log.in <- read.csv("Ent_2013_logger_data.csv")

  setwd(paste0(mainPath,yearPath))

  Log.in$SiteName <- as.character(Log.in$SiteName)
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)

  LST.Log.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteName  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$Entiat_RCA[ID.in$SiteName == i]
      Elev <- ID.in$Elev_M[ID.in$SiteName == i]
      
      LST.site <- matrix(ncol=3, nrow=46)
      LST.site[,1] <- newnamesnum
      LST.site[,2] <- unlist(LST.in[LST.in$RCAID == RCAID,1:46])
      LST.site <- data.frame(LST.site)
      colnames(LST.site) <- c("yrJul", "LST", "Elev")
      LST.site[3] <- Elev
      LST.Log.site <- merge(LST.site, Log.site, by.x = "yrJul", by.y = "yrJul", all.x=TRUE, all.y = FALSE)
      
      LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
    }


  ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
  NoNA.xyz <- LST.Log.out[ind,]

  NoNA.xyz <- NoNA.xyz[,c(18, 2, 1, 3, 7)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  plot(NoNA.xyz$x, NoNA.xyz$y)

  write.csv(x=NoNA.xyz, file=paste0("Ent_2013_8Day_Min_model_data.csv", row.names = FALSE)

  NoNA.xyz <- orderBy(~z, NoNA.xyz)

  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:maxrow,]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]

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
  MSS <- sum((pred.y - mean(y))^2)
  p2 <- MSS/TSS

  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "summer"
  pred.out[,5] <- "2013"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2013_summer_only.csv",sep = ",", col.names=F)  
  
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
  MSS <- sum((pred.y - mean(y))^2)
  p2 <- MSS/TSS

  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "full year"
  pred.out[,5] <- "2013"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2013_full_year.csv",sep = ",", col.names=F)  
  
  plot(pred.out[,1], pred.out[,2])
  summer_pred <- subset(pred.out, z > 181 & z < 258)
  points(summer_pred[, 1], summer_pred[,2], pch = 16, col = "green")
  abline(0,1)
         
#####################################
# spring/fall
######################################

coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bJul=numeric(2), bElev=numeric(2))
metrics_out <- data.frame(PRESS=numeric(2), p2=numeric(2), r2=numeric(2), RMSEP=numeric(2), RMSE=numeric(2))
rownames(metrics_out) <- c("Spring", "Fall")
rownames(coeffs_out) <- c("Spring", "Fall")


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
  RSS <- sum(pressstat_sum$residuals^2)
  TSS <- sum((y - mean(y))^2)
  MSS <- sum((pred.y - mean(y))^2)
  p2 <- MSS/TSS

  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "spring"
  pred.out[,5] <- "2013"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_End_2013_sp_fall.csv",sep = ",", col.names=F)  

  sp_pred <- subset(pred.out, z > 181 & z < 258)
  points(sp_pred[, 1], sp_pred[,2], pch = 16, col = "red")


  coeffs_out[1,1] <- coeffs[1,1]
  coeffs_out[1,2] <- coeffs[2,1]
  coeffs_out[1,3] <- coeffs[3,1]
  coeffs_out[1,4] <- coeffs[4,1]
  coeffs_out[1,5] <- coeffs[5,1]
  
  metrics_out[1,1] <- pressstat_sum$stat
  metrics_out[1,2] <- p2
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



  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
  RSS <- sum(pressstat_sum$residuals^2)
  TSS <- sum((y - mean(y))^2)
  MSS <- sum((pred.y - mean(y))^2)
  p2 <- MSS/TSS

  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "fall"
  pred.out[,5] <- "2013"

  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2013_sp_fall.csv",sep = ",", col.names=F)  

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

write.table(x=coeffs_out, append=F,row.names=T, file = "All_data_2013_mod_coeffs_Mx.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = "All_data_2013_mod_metrics_Mx.csv", sep = ",", col.names=T)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################
setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Entiat/")
elev.in <- read.csv("Ent_rca_elev.csv")
LST.in <- read.dbf("LST13_Ent_RCA.dbf")

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
yearPath <- "2013"

  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- substring (colnames(LST.in)[1:46],3,5)
  newnamesnum <- as.numeric(newnamesnum)
  colnames(LST.in)[1:46] <- newnamesnum
  
  LST.elev <- merge(LST.in, elev.in, by.x = "Ent_RCAID", by.y = "rca_id")
  
  setwd(paste0(mainPath,yearPath))
  
  coeffs.in <- read.csv("All_data_2013_mod_coeffs_Mx.csv")
  LogPred.out <- LST.elev[,c(25:36,48,1)]
  LogPred.out[,1:12] <- 0
  rcas <- unique(elev.in$rca_id)
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
LogPred.out$Basin_RCA <- paste0("Ent_", LogPred.out$Ent_RCAID)
namesnum <- as.numeric(colnames(LogPred.out[1:12]))

names.out <- sprintf("Tmx_13_%03d", namesnum)
colnames(LogPred.out)[1:12] <- names.out[1:12]


write.dbf(LogPred.out, file = "predt2013_Ent_8D_Max_summer.dbf") 


#______________________
# went back and calc-ed the model metrics for sping & fall together:

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/2013/")

Log.in <- read.csv("Ent_2013_8Day_model_data.csv")
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

pred.new <- predict(mod, data = data.sp)
pred.new[pred.new < -0.5] = -0.5
data.sp$pred <- unlist(pred.new)

pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")

resids <- pressstat_sum$residuals

y <- data.fall$y
x <- data.fall$x
z <- data.fall$z
e <- data.fall$e
mod <- lm(y ~ x + I(x^2) + z + e)
sum_mod <- summary(mod)

pred.new <- predict(mod, data=data.fall)
pred.new[pred.new < -0.5] = -0.5
data.fall$pred <- unlist(pred.new)

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
# summary for the summer max EPs linear fill
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
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

  Au.in <- read.csv("Ent_RCA_AUs_IP.csv")
  ST_rcas <- Au.in[Au.in$STHD_IP == "Steelhead",]
  CH_rcas <- Au.in[Au.in$CHNK_IP == "Chinook",]

  setwd(paste0(mainPath,yearPath))

  Max.in <- read.dbf("predt2013_Ent_8D_Max_summer.dbf") 

#____some shenanigans to fill in daily temps so averages over time window works_____

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

  filled.out$RCAID <- Max.in$Ent_RCAID
  colnames(filled.out)[1:49] <- 201:249
#__________________________________________________

  Max.au <- merge(filled.out, Au.in, by.x = "RCAID", by.y = "rca_id")

  Max.au$STHD_IP <- as.character(Max.au$STHD_IP)
  Max.au$CHNK_IP <- as.character(Max.au$CHNK_IP)
  Max.au$AU_ID <- as.character(Max.au$AU_ID)

  


#_________all_______________________  
  rcas <- unique(Max.au$RCAID)
  
  SumSumm.out <- data.frame ("RCAID" = rcas)


# count of days in exceedence -----------------------------------------


  for (i in 1:length(rcas)) 
    { 
      l <- rcas[i]
      MaxRCA <- Max.au[Max.au$RCAID == l,] #grab days for one RCA 
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


SumSumm.out[,2:7] <- round(SumSumm.out[,2:7], digits = 2)

colnames(SumSumm.out) <- c("RCAID", "Pct12_2013","Pct13_2013", "Pct16_2013", "Pct18_2013", "Pct20_2013", "Pct22_2013",  "MxMx_2013", "sdMn_2013", "MnMx_2013")

write.dbf(SumSumm.out, file = "Ent_2013_21Jul_31Aug_max_summary_All.dbf")
write.csv(SumSumm.out, file = "Ent_2013_21Jul_31Aug_max_summary_All.csv", row.names = F)


#________________________________
# summmary of IP reaches only by AU
#______________________________________


SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "rca_id")
colnames(SumSummAU)[1:10] <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx")
SumSummAU$STHD_IP <- as.character(SumSummAU$STHD_IP)
SumSummAU$CHNK_IP <- as.character(SumSummAU$CHNK_IP)
SumSummAU$AU_ID <- as.character(SumSummAU$AU_ID)

#________steelhead____________________  

Aus <- unique(Max.au$AU_ID)
rcas <- unique(ST_rcas$rca_id)



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
      SumAU <- SumSummAU[SumSummAU$AU_ID == Aus[i] & SumSummAU$STHD_IP == "Steelhead",]
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
write.csv(SumAU.out, file = "Ent_AU_2013_21Jul_31Aug_ExPcnt_summary_STHD_IP.csv", row.names = F)


#________chinook_______________________
Aus <- unique(Max.au$AU_ID)
rcas <- unique(CH_rcas$rca_id)

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


  for (i in 1:length(Aus)) 
    { 
      SumAU <- SumSummAU[SumSummAU$AU_ID == Aus[i] & SumSummAU$CHNK_IP == "Chinook",]
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
write.csv(SumAU.out, file = "Ent_AU_2013_21Jul_31Aug_ExPcnt_summary_CHNK_IP.csv", row.names = F)

# #####_________________pie charts_____________


library(plotrix)

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
  yearPath <- "2013"
  setwd(paste0(mainPath, yearPath))
  
  metric <- "22"
  popn <- "CHNK"
  varName <- paste0("v", metric, "DMax")

  SumAU.out <- read.csv(paste0("Ent_AU_2013_21Jul_31Aug_ExPcnt_summary_", popn, "_IP.csv"))
  SumAU.out$AU_Code <- as.character(SumAU.out$AU_Code)
  Aus <- SumAU.out$AU_Code
  Aus <- Aus[-5]

 
  for (i in 1:length(Aus)) 
    { 
      ExVal <- SumAU.out[SumAU.out$AU_Code == Aus[i], varName]*100
      value <- c(ExVal, 100-ExVal)
      name <- as.character(Aus[i])
      
      filename <- paste0("D:/OneDrive/work/research/CHaMP/graphics/EP/Ent/2013/PctD", metric, "/", name, "_", metric, popn, ".png", sep="")
      png(filename=filename)
      
      
      if(ExVal >= 10)
      {
        cols <- c("red2", "gainsboro")
        pie(value, col=cols, cex.main=2.0,)
        bisect.angles <- floating.pie(0,0,value, col=cols)
        pie.labels(0,0,bisect.angles,radius=0.4, c(paste0(value[1],"%")), cex=3, font=2, main=name)
        
      } else if (ExVal < 10 & ExVal > 0){
        cols <- c("red2", "gainsboro")
        pie(value, col=cols, cex.main=2.0,)
        bisect.angles <- floating.pie(0,0,value, col=cols)
        pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=3, font=2, main=name)
      } else {
        cols <- c("red2", "gainsboro")
        pie(value, clockwise = TRUE, col=cols, cex.main=2.0,)
        bisect.angles <- floating.pie(0,0,value, col=c("gainsboro"), )
        pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=3, font=2, main=name)
      }
      
      dev.off()
    }

# #########################
#This parts formats the error by day/site info
# ##########################
  
  error.pts <- rbind(data.sp, data.fall)
  error.pts$SiteName <- as.character(error.pts$SiteName)
  SiteID <- unique(error.pts$SiteName)
  SiteID <- as.matrix(SiteID)
  Error.pts.out <- matrix(nrow = length(SiteID), ncol = (length(names.out)+1))
  
  errorName <- paste0("JulDay_", namesnum)
  
  Error.pts.out[1:length(SiteID), 1] <- unlist(SiteID)[1:length(SiteID)]
  Error.pts.out <- as.data.frame(Error.pts.out, stringsAsFactors = FALSE)
  colnames(Error.pts.out)[2:(length(names.out)+1)] <- errorName
  colnames(Error.pts.out)[1] <- "SiteName"
  
  for (i in SiteID) 
    { 
      error.site <- error.pts[error.pts$SiteName  == i,]
      error.site$error <- error.site[,'y']-error.site[,'pred']
      
      
      error <- matrix(ncol=1, nrow=length(names.out))
      error[,1] <- namesnum
      error <- data.frame(error)
      colnames(error) <- c("JulDay")
      
      error.site.fill <- merge(error, error.site, by.x = "JulDay", by.y = "z", all.x=TRUE, all.y = FALSE)
      Error.pts.out[Error.pts.out$SiteName==i,2:(length(names.out)+1)] <- as.numeric(unlist(error.site.fill$error))
    }
  
  Error.pts.out[,2:(length(names.out)+1)] <- sapply(Error.pts.out[,2:(length(names.out)+1)], as.numeric)
  
  Error.pts.out[,2:(length(names.out)+1)] <- round(Error.pts.out[,2:(length(names.out)+1)], digits=3)
  
  write.dbf(Error.pts.out, file = paste0("Error", "_", basin, "_", yearPath, "_8D_Mn.dbf")) 
  write.csv(Error.pts.out, file = paste0("Error", "_", basin, "_", yearPath, "_8D_Mn.csv"))
  
# ##############################################
# shapefile output
# ##############################################
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
  
  basin <- "Ent"
  midBasin <- "Entiat"
  longBasin <- "Entiat"
  yrPath <- "13"
  yearPath <- "2013"
  subDir <- "Entiat/"
  dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)
  
  projname <- paste0(longBasin, "_net_rca")
  ptname <- paste0(basin, "_sites_rca_pj")
  setwd(paste0(dataPath, longBasin))
  
  proj_layer <- readOGR(dsn=".", layer = projname)
  
  netname <- paste0(basin, "_STHD_net")
  network <- readOGR(dsn=".", layer = netname)
  network <- spTransform(network, proj4string(proj_layer))
  
  setwd(paste0(mainPath, longBasin, "/", yearPath, "/"))
  
  preds <- read.dbf(paste0("predt", yearPath, "_", basin, "_8D_Mn.dbf"))
  netmerge <- merge(network, preds, by.x='RCAID', by.y = 'RCAID')
  
  
# #### plot it to make sure it looks right ######### 
  
  
  fix3 <- classIntervals(netmerge@data[,10], n = 11, style = "fixed",fixedBreaks=c(-1,0,1,2,4,6,8,10,12,14,20))
  fix3.colors <- findColours(fix3,pal=seis)
  plot(netmerge, col=fix3.colors, bg="black", fg="white")
  
# #################################################  
  
  setwd(paste0(dataPath, longBasin))
  layername <- paste0(basin, "_", yearPath, "_8D_mn")
  writeOGR(netmerge, dsn=".", layer=layername, driver="ESRI Shapefile")
  
###############################################
#
# Edited 5 may 2016 to add animation output
################################################
  
  library(rgdal)
  library(RColorBrewer)
  library(classInt)
  
  shpPath <- paste0(modelPath, "graphics/")
  modelPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", yearPath, "_temp_CHaMP/")
  ptsPath <- paste0("D:/OneDrive/work/research/CHaMP/CHaMP_data/", longBasin, "/")
  
  netname <- paste0(basin, "_", yearPath, "_8D_mn")
  ptsname <- paste0(basin, "_Error_", yearPath, "_8D_mn")
  setwd(ptsPath)
  
  network <- readOGR(dsn=paste0(mainPath, yearPath, "_temp_CHaMP"), layer = netname)
  error_pts <- readOGR(dsn=".", "Ent_Error_2014_8D_Mn")
  
  setwd(shpPath)
  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)
  
  
  names.out <- colnames(network@data[3:48])
  namesnum <- as.numeric(gsub("Tmn_14_", "", colnames(network@data[3:48])))
  means <- colMeans(network@data[3:48])
  SDs <- colStdevs(network@data[3:48])
  yplus <- means + SDs
  yminus <- means - SDs
  df <- data.frame(means=means, SDs=SDs, names=namesnum)
  sequ <- c(1:46)
  namer <- sprintf('%03d', sequ)
  fix4 <- classIntervals(means, n = 10, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18))
  fix4.colors <- findColours(fix4,pal=seis)
  
  for (i in 3:48)
  {
    
    namey <- gsub("Tmn_14_", "", colnames(network@data)[i])
    
    filename <- paste0(mainPath, longBasin, "/", yearPath, "/graphics2/", namer[i-3], ".png", sep="")
    png(filename=filename, res = 300, width = 1500, height = 1500, units = "px", bg="black")
    
    fix3 <- classIntervals(network@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,2,4,6,8,10,12,14,16,18,22))
    fix3.colors <- findColours(fix3,pal=seis)
    
    cexEr <- ifelse(abs(error_pts@data[,i-1]) <= 1, 0.5,
                    ifelse(abs(error_pts@data[,i-1])>1, 1,
                           ifelse(abs(error_pts@data[,i-1])>2, 1.5,
                                  ifelse(abs(error_pts@data[,i-1])>3, 2.0, NA))))
    
    plot(network, col=fix3.colors, bg="black", fg="white")
    points(error_pts, pch=21, col="black", bg="gray40", cex=cexEr)
    
    legend("right", fill = attr(fix3.colors, "palette"), title="8-day Mn (°C)", legend = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","20-22"), bty = "n", cex=.5, inset=c(.1,0), text.col="white");
    legend(x=grconvertX(c(0.05, 0.25), from='npc'), 
           y=grconvertY(c(0.3, 0.5), from='npc'),, title="Model error (°C)", legend = c("0-1","1-2","2-3","3+"), bty = "n", pch=16, pt.cex=c(0.5, 0.75, 1.0, 1.5), col = "gray40", cex=.5, text.col="white");
    
    
    title("Entiat 8-day mean 2014 (°C)", sub = paste0("Julian Day ", namey), line=-0.9, adj=.80, col.main="white", col.sub="white", outer=FALSE, cex.main=0.5, cex.sub=0.5)
    tmp2 <- subplot(
      plot(namesnum[1:(i-3)], means[1:(i-3)], col=fix4.colors, pch=16, bty="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean (+/-SD)", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"), 
      x=grconvertX(c(0.1,0.45), from='npc'), 
      y=grconvertY(c(0.05, 0.20), from='npc'),
      size=c(1,1.5), vadj=0.5, hadj=0.5, 
      pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
    op <- par(no.readonly=TRUE)
    par(tmp2)
    arrows(namesnum[1:(i-3)], yplus[1:(i-3)], namesnum[1:(i-3)], yminus[1:(i-3)], length=0, lwd=5, code=3, lend=0, col="gray20")
    par(op)
    tmp2 <- subplot(
      plot(namesnum[1:(i-3)], means[1:(i-3)], col=fix4.colors, pch=16, bty="n", xlim=c(0,360), ylim=c(0,18), cex.main=.8, main="Basin mean (+/-SD)", adj=0, xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"), 
      x=grconvertX(c(0.1,0.45), from='npc'), 
      y=grconvertY(c(0.05, 0.20), from='npc'),
      size=c(1,1.5), vadj=0.5, hadj=0.5, 
      pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
    
    dev.off()
  }
  
  
  setwd(paste0(mainPath, longBasin, "/", yearPath, "/graphics2/"))
  
  system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png Entiat_2014_8day_mn.mpeg')
  
  