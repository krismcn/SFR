############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for a summer max temps MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# 
# Edited Aug 2014 to add the PRESS stastic output
# Should look at PRESS output to screen to chose model - only that model will output resids and preds.

#Edited Oct 2015 for EP watershed summer output
#Edited Jan 2016 to update the gap-filling interpolation functions
          
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
library(zoo)

  basin <- "CW"
  medBasin <- "Cwater"
  longBasin <- "Clearwater"
  year <- "11"
  yearPath <- "2011"

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/"
  
  setwd(paste0(mainPath, longBasin))
           
####################################
# Parsing out the subwatersheds
######################################

  NoNA.xyz <- read.csv(paste0(basin, "_", yearPath, "_8Day_model_data.csv"), stringsAsFactors=FALSE)
  
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  
  
  basin <- "LCW"
  setwd(paste0(mainPath, longBasin))

  Sub_ID.in <- read.dbf(paste0(basin,"_sites_elev.dbf"), as.is=TRUE)
  
  Sub_IDs <- unique(Sub_ID.in$Orig_ID)

  #Sub_data <- summer[summer$SiteName %in% Sub_IDs, ]
  #Sub_data$y <- log(Sub_data$y)



  Sub_data <- NoNA.xyz[NoNA.xyz$SiteName %in% Sub_IDs, ]
  #Sub_data$y[Sub_data$y < 0.1] = 0.1
  #Sub_data$y <- log(Sub_data$y)

  y <- Sub_data$y
  x <- Sub_data$x
  z <- Sub_data$z
  e <- Sub_data$e
  plot(z, y)  
  plot(x, y)
  
  mod <- lm(y ~ x + I(x^2) + e)
  sum_mod <- summary(mod)
  pred.y <- predict(mod)
  plot(y, pred.y)

  maxrow <- which.max(Sub_data$y)
  data.sp <- Sub_data[1:maxrow,]
  data.fall <- Sub_data[maxrow:nrow(Sub_data),]

  
  
  
#################################
# to suss out summer stuff
###############################

summer <- subset(Sub_data, z > 181 & z < 258)
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

setwd(paste0(mainPath, longBasin, "/", yearPath))
write.table (x=pred.out,append=T,row.names=F,file= paste0("jk_pred_v_y_", basin, "_", yearPath, "_summer_only.csv"),sep = ",", col.names=F)  

plot(pred.out[,1], pred.out[,2])
abline(0,1)

################################
# full year
###################################

  setwd(paste0(mainPath, longBasin, "/", yearPath))

#   y <- NoNA.xyz$y
#   x <- NoNA.xyz$x
#   z <- NoNA.xyz$z
#   e <- NoNA.xyz$e
#   plot(x, y)
#   


  mod <- lm(y ~ x + I(x^2) + e)
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
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_full_year.csv"),sep = ",", col.names=F)  
  
  plot(pred.out[,1], pred.out[,2])
  summer_pred <- subset(pred.out, z > 181 & z < 258)
  points(summer_pred[, 1], summer_pred[,2], pch = 16, col = "green")
  abline(0,1)

#####################################
# spring/fall
######################################
  setwd(paste0(mainPath, longBasin, "/", yearPath))

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
  pred.out[,5] <- yearPath
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall.csv"), sep = ",", col.names=F)  

  sp_pred <- subset(pred.out, z > 181 & z < 258)
  points(sp_pred[, 1], sp_pred[,2], pch = 16, col = "red")
  plot(sp_pred[, 1], sp_pred[,2], pch = 16, col = "red")


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
  pred.out[,5] <- yearPath

  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  

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

write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mx.csv"), sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_Mx.csv"), sep = ",", col.names=T)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################
  
  basin <- "CW"
  setwd(paste0(mainPath, longBasin))
  elev.in <- read.csv(paste0(basin, "_rca_elev.csv"))
  LST.in <- read.dbf(paste0("LST", year, "_", medBasin, "_RCA.dbf"), as.is=TRUE)


  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));
  colnames(LST.in)[1:46] <- newnamesnum;

  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")


  basin <- "LCW"
  setwd(paste0(mainPath, longBasin, "/", yearPath))

  coeffs.in <- read.csv(paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mx.csv"))
  LogPred.out <- LST.elev[,c(25:36,48,1)]
  LogPred.out[,1:12] <- 0
  rcas <- unique(elev.in$rca_id)
  LST.sum <- LST.elev[,c(25:36, 48,1)]

  for (i in 1:length(rcas))  
    {
      x <- unlist(LST.sum[i,])
      maxrow <- as.numeric(which.max(x[1:12])) #either specify or let be dynamic
      midrow <- maxrow - 1
      day <- as.numeric(colnames(LST.sum)[maxrow])
      
      j <- 185
      for (l in 1:maxrow)
      {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
       j <- j + 8}
      k <- day
      for (l in midrow:12)     
      {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[13] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
       k <- k + 8}
      LogPred.out[i,1:12] <- x [1:12] 
    }

#   for (i in 1:length(rcas))  
#     {
#       x <- unlist(LST.sum[i,])
#       
#       j <- 185
#       for (l in 1:12)
#       {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
#        j <- j + 8}
#       
#       LogPred.out[i,1:12] <- x [1:12] 
#     }

  LogPred.out <- as.data.frame(LogPred.out)

  #LogPred.out[,1:12] <- exp(LogPred.out[,1:12])

  LogPred.out[LogPred.out< -0.5] = -0.5
  LogPred.out$Basin_RCA <- paste0(basin, "_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:12]))
  varName <- paste0("Tmx_", year)
  names.out <- sprintf("%s_%03d", varName, namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]


write.dbf(LogPred.out, file = paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf")) 



#________________________________________________________
# summary for the summer max EPs linear fill

############################
# 15Jul-31Aug ------------------------------------------------------------------
#############################

  setwd(paste0(mainPath, longBasin))
  basin <- "CW"
  popn <- "STHD"
  varName <- paste0(basin, "_", popn)

  Au.in <- read.csv(paste0(basin, "_RCA_AUs_IP.csv"), stringsAsFactors = FALSE)
  
  Aus <- unique(Max.au$LCW_AU)
  Aus <- Aus[-3]
  rcas <- unique(ST_rcas$rca_id)

  basin <- "LCW"

  blerg <- paste0(basin, "_STHD")
  blarg <- paste0('ST_rcas <- Au.in[Au.in$', blerg, ' == "Steelhead",]')
  eval(parse(text = blarg))

  
  setwd(paste0(mainPath, longBasin, "/", yearPath))

  blee <- paste0("Tmx_", year, "_")
  bloo <- paste0('colnames(Max.in) <- gsub("', blee, '", "", colnames(Max.in))')
  

  Max.in <- read.dbf(paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf"))
  eval(parse(text = bloo))

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




#_________all_______________________  
rcas <- unique(Max.au$RCAID)

SumSumm.out <- data.frame ("RCAID" = rcas)


# count of days in exceedence -----------------------------------------


  for (i in 1:length(rcas)) 
    { 
      l <- rcas[i]
      MaxRCA <- Max.au[Max.au$RCAID == l,] #grab days for one RCA 
      DaysAbove12 <- length(which(MaxRCA[2:44]> 12)) #finds how many day in the 20July-31Aug window exceed threshold 
      DaysAbove13 <- length(which(MaxRCA[2:44]> 13))
      DaysAbove16 <- length(which(MaxRCA[2:44]> 16))
      DaysAbove18 <- length(which(MaxRCA[2:44]> 18))
      DaysAbove20 <- length(which(MaxRCA[2:44]> 20))
      DaysAbove22 <- length(which(MaxRCA[2:44]> 22))
      MaxMax <- max(MaxRCA[2:44])
      SDMax <- sd(MaxRCA[2:44])
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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2011","Pct13_2011", "Pct16_2011", "Pct18_2011", "Pct20_2011", "Pct22_2011",  "MxMx_2011", "sdMn_2011", "MnMx_2011")

setwd(paste0(mainPath, longBasin, "/", yearPath))
write.dbf(SumSumm.out, file = paste0(basin,"_", yearPath, "_21Jul_31Aug_max_summary_All.dbf"))
write.csv(SumSumm.out, file = paste0(basin,"_", yearPath, "_21Jul_31Aug_max_summary_All.csv"), row.names = F)


#________________________________
# summmary of IP reaches only by AU
#______________________________________


SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "rca_id")
colnames(SumSummAU)[1:10] <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx")


######
#________steelhead____________________  
######

  blerg <- paste0(basin, "_AU")
  blarg <- paste0('Aus <- unique(Max.au$', blerg, ')')
  eval(parse(text = blarg))


#######_____Need to generalize this section_________#####

  Aus <- unique(Max.au$LCW_AU)
  Aus <- Aus[-3]
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


#######_____Need to generalize this section_________#####

  for (i in 1:length(Aus)) 
    { 
      SumAU <- SumSummAU[SumSummAU$LCW_AU == Aus[i] & SumSummAU$LCW_STHD == "Steelhead",]
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
write.csv(SumAU.out, file = paste0(basin, "AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, "_IP.csv"), row.names = F)




#######
#_________________pie charts_____________
#######
`

library(plotrix)

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/"
yearPath <- "2011"
setwd(paste0(mainPath, longBasin, "/", yearPath))
graphPath <- "D:/OneDrive/work/research/CHaMP/graphics/EP/"
basin <- "LCW"
popn <- "STHD"


  SumAU.out <- read.csv(paste0(basin, "AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, "_IP.csv"), stringsAsFactors = FALSE)
  
  Aus <- SumAU.out$AU_Code


  metric <- "18"
  varName <- paste0("v", metric, "DMax")

  for (i in 1:length(Aus)) 
    { 
      ExVal <- SumAU.out[SumAU.out$AU_Code == Aus[i], varName]*100
      value <- c(ExVal, 100-ExVal)
      name <- as.character(Aus[i])
      
      filename <- paste0(graphPath, basin, "/", yearPath, "/PctD", metric, "/", name, "_", metric, popn, ".png", sep="")
      png(filename=filename)
      
      
      if(ExVal >= 10)
      {
        cols <- c("red2", "gainsboro")
        pie(value, col=cols, cex.main=3.0, main=name)
        bisect.angles <- floating.pie(0,0,value, col=cols,)
        pie.labels(0,0,bisect.angles,radius=0.4, c(paste0(value[1],"%")), cex=3, font=2, main=name)
        
      } else if (ExVal < 10 & ExVal > 0){
        cols <- c("red2", "gainsboro")
        pie(value, col=cols, cex.main=3.0,main=name)
        bisect.angles <- floating.pie(0,0,value, col=cols,)
        pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=3, font=2, main=name)
      } else {
        cols <- c("red2", "gainsboro")
        pie(value, clockwise = TRUE, col=cols, cex.main=3.0,main=name)
        bisect.angles <- floating.pie(0,0,value, col=c("gainsboro"), )
        pie.labels(0,0,bisect.angles,radius=1.05, c(paste0(value[1],"%")), cex=3, font=2, main=name)
      }
      
      dev.off()
    }