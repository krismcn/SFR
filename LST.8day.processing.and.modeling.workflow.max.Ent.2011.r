############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for a year using MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# Currently set up for the John Day which has 12039 grid cells and 5532 RCAs (those parameters can be changed for other regions).
# Should do a search-and-replace for the output folder (usually dated), and the year being processed, in both a YYYY and a _YY format.

# Edited Aug 2014 to add the PRESS stastic output
# Should look at PRESS output to screen to chose model - only that model will output resids and preds.


          
##############################################################################################
# This section reads in the 1km LST data for a year, uses a 4th order polynomial to fill in Julian day 1 & 365 (if they are missing)
# then fills any remaining gaps across the year at each pixel with a linear interpolation.
##############################################################################################
setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Entiat")

library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)
library(zoo)

  LST.in <- read.csv("Ent_LST_2011.csv", header = TRUE)
  GrPolID <- LST.in[,1]
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
  LST.out$GRID_CODE <- GrPolID
  
  newnames <- substring (colnames(LST.out)[1:46],11,13); #clips the MODIS grid names down to the julian date
  newnamesnum <- as.numeric(newnames);    #specifies those dates to be numeric in order to drop the leading zeros
  colnames(LST.out)[1:46] <- newnamesnum;
  
  write.dbf(LST.out, file = "LST12_Ent_interp.dbf")

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

out_Preds <- read.dbf("LST11_Ent_interp.dbf") #use the appropriate read statement

weights <- read.dbf("Ent_rca_area_wgts.dbf")
weights[weights$area_wgts<0.01, "area_wgts"] = 0.01
rca_zonal <- matrix(ncol=46, nrow=448) #change for the number of actual RCAs in the watershed

for(i in 1:448)	#assumes the RCAs are consecutively numbered starting with "1"
    {
    pixels <- weights[weights$rca_id == i, "GRIDCODE_1"]
    wgts <- weights[weights$rca_id == i, "area_wgt"]
    
    for (j in 1:46)
        {
        daily <- out_Preds[out_Preds$GRID_CODE %in% pixels, j]
        zonal_mn <- weighted.mean(daily, wgts, na.rm = TRUE)
        rca_zonal[i,j] <- zonal_mn
        }
    }


colnames(rca_zonal)[1:46] <- colnames(out_Preds)[1:46]
rca_zonal <- as.data.frame(rca_zonal)
rca_zonal$RCAID <- 1:448

write.dbf(rca_zonal, file = "LST11_Ent_RCA.dbf")



#################################################################
#Logger prediction modeling part

#################################################################
setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Entiat")

  ID.in <- read.csv("Ent_sites_elev.csv")


  LST.in <- read.dbf("LST11_Ent_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(colnames(LST.in)[1:46]);    
  colnames(LST.in)[1:46] <- newnamesnum;


mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
yearPath <- "2011"
setwd(paste0(mainPath))

  Log.in <- read.csv("Ent_2011_logger_data.csv")

  setwd(paste0(mainPath,yearPath))

  Log.in$SiteName <- as.character(Log.in$SiteName)
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)
  ID.in$SiteName <- as.character(ID.in$SiteName)

  LST.Log.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteName  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$Entiat_RCA[ID.in$SiteName == i]
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

  NoNA.xyz <- NoNA.xyz[,c(15, 2, 1, 3, 7)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  plot(NoNA.xyz$x, NoNA.xyz$y)

  write.csv(x=NoNA.xyz, file="Ent_2011_8Day_model_data.csv", row.names = FALSE)

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
  pred.out[,5] <- "2011"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2011_summer_only.csv",sep = ",", col.names=F)  
  
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
  pred.out[,5] <- "2011"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2011_full_year.csv",sep = ",", col.names=F)  
  
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
  pred.out[,5] <- "2011"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_End_2011_sp_fall.csv",sep = ",", col.names=F)  

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
  pred.out[,5] <- "2011"

  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2011_sp_fall.csv",sep = ",", col.names=F)  

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

write.table(x=coeffs_out, append=F,row.names=T, file = "All_data_2011_mod_coeffs_Mx.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = "All_data_2011_mod_metrics_Mx.csv", sep = ",", col.names=T)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################
setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Entiat/")
elev.in <- read.csv("Ent_rca_elev.csv")
LST.in <- read.dbf("LST11_Ent_RCA.dbf")

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
yearPath <- "2011"

colnames(LST.in)<-gsub("X", "", colnames(LST.in))
newnamesnum <- as.numeric(colnames(LST.in)[1:46])
colnames(LST.in)[1:46] <- newnamesnum;

LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")

setwd(paste0(mainPath,yearPath))

coeffs.in <- read.csv("All_data_2011_mod_coeffs_Mx.csv")
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
LogPred.out$Basin_RCA <- paste0("Ent_", LogPred.out$RCAID)
namesnum <- as.numeric(colnames(LogPred.out[1:12]))

names.out <- sprintf("Tmx_11_%03d", namesnum)
colnames(LogPred.out)[1:12] <- names.out[1:12]


write.dbf(LogPred.out, file = "predt2011_Ent_8D_Max_summer.dbf") 


#______________________
# went back and calc-ed the model metrics for sping & fall together:

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/2011/")

Log.in <- read.csv("Ent_2011_8Day_model_data.csv")
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
# summary for the summer max EPs linear fill
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
yearPath <- "2011"


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
  
  Max.in <- read.dbf("predt2011_Ent_8D_Max_summer.dbf") 
  
  colnames(Max.in)<-gsub("Tmx_11_", "", colnames(Max.in))

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

  Max.au <- merge(Max.out, Au.in, by.x = "RCAID", by.y = "rca_id")

#____________________________  
  Aus <- unique(Max.au$AU_STHD)
  rcas <- unique(ST_rcas$RCA_ID)
#_______________________________
  Aus <- unique(Max.au$AU_CHNK)
  rcas <- unique(CH_rcas$RCA_ID)
#________________________________  
  rcas <- unique(Max.au$RCAID)

SumSumm.out <- data.frame ("RCAID" = rcas)



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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2011","Pct13_2011", "Pct16_2011", "Pct18_2011", "Pct20_2011", "Pct22_2011",  "MxMx_2011", "sdMn_2011", "MnMx_2011")

write.dbf(SumSumm.out, file = "Ent_2011_21Jul_31Aug_max_summary_All.dbf")
write.csv(SumSumm.out, file = "Ent_2011_21Jul_31Aug_max_summary_All.csv", row.names = F)


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
write.csv(SumAU.out, file = "Ent_AU_2011_21Jul_31Aug_ExPcnt_summary_STHD_IP.csv", row.names = F)


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
write.csv(SumAU.out, file = "Ent_AU_2011_21Jul_31Aug_ExPcnt_summary_CHNK_IP.csv", row.names = F)

#_________________pie charts_____________


library(plotrix)

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
  yearPath <- "2011"
  setwd(paste0(mainPath, yearPath))
  
  popn <- "CHNK"
  metric <- "18"
  varName <- paste0("v", metric, "DMax")

  
  SumAU.out <- read.csv(paste0("Ent_AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, "_IP.csv"))
  SumAU.out$AU_Code <- as.character(SumAU.out$AU_Code)
  Aus <- SumAU.out$AU_Code
  Aus <- Aus[-5]

  
    for (i in 1:length(Aus)) 
      { 
        ExVal <- SumAU.out[SumAU.out$AU_Code == Aus[i], varName]*100
        value <- c(ExVal, 100-ExVal)
        name <- as.character(Aus[i])
        
        filename <- paste0("D:/OneDrive/work/research/CHaMP/graphics/EP/Ent/", yearPath, "/PctD", metric, "/", name, "_", metric, popn, ".png", sep="")
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
