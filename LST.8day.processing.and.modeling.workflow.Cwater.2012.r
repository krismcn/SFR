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
  year <- "12"
  yearPath <- "2012"

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/"
  setwd(paste0(mainPath, "shapes/"))

 # LST.in <- read.dbf("CW_11_LST.dbf")
  #write.csv(x = LST.in, file = "CW_11_LST.csv", row.names = F)
  #LST.in <- LST.in[,1:2]
  #write.dbf(LST.in, "CW_11_LST.dbf")
  #fixed <- read.csv("CW_11_LST.csv")
  #write.dbf(fixed, "CW_11_LST.dbf")  

  #LST.in <- read.dbf(paste0(basin, "_", year, "_LST.dbf"))
  LST.in <- read.csv(paste0(basin, "_", year, "_LST.csv"))

  PointID <- LST.in[,1]
  

#-------------------this part fills in some holes along the seam between swath 1 & 2---------
  LST.in[LST.in$POINTID == 608, 3:48] <- LST.in[LST.in$POINTID == 609, 3:48]
  LST.in[LST.in$POINTID == 656, 3:48] <- LST.in[LST.in$POINTID == 657, 3:48]
  LST.in[LST.in$POINTID == 704, 3:48] <- LST.in[LST.in$POINTID == 705, 3:48]
  LST.in[LST.in$POINTID == 14990, 3:48] <- LST.in[LST.in$POINTID == 9276, 3:48]
  LST.in[LST.in$POINTID == 14991, 3:48] <- LST.in[LST.in$POINTID == 9557, 3:48]
  LST.in[LST.in$POINTID == 14992, 3:48] <- LST.in[LST.in$POINTID == 9838, 3:48]
  LST.in[LST.in$POINTID == 14993, 3:48] <- LST.in[LST.in$POINTID == 10116, 3:48]

  LST.in <- LST.in[,3:48]
  LST.in[LST.in<0.1] = NA

##############################
#_____creates a data frame using a 3 window rollaply to fill in rows with too many NAs  
 
# tLST.in <- t(LST.in)
# tLST.out <- rollapply(tLST.in, width = 4, FUN = mean, na.rm = FALSE, align = "left")

# tCLST <- t(tLST.out)
# CLST.out <- as.data.frame(tCLST)
# CLST <-c(CLST.out[,1:43], LST.in[,44:46])
# CLST <- as.data.frame(CLST)

#______to use the 3-day rolling mean as gap filling_______________



#__________This function doesn't work, but I'm leaving it in here as a record______________

#   myFunc <- function(a,b)
#     {
#       x <- unlist(a)
#       y <- unlist(b)
#       
#       for (l in 1:46)
#       {
#         if(is.na(x[l]))
#         {
#           x[l] <- y[l] 
#         }   
#       }
#      a <- x
#      return(a)
#     }
# 
#   LST.filled2 <- apply(LST.test, 1, function(x) myFunc(LST.test, CLST))
#  
#___This function worked, but took a long time to process (30+ minutes)_____
#___Leaving it in for a record________________________________
#   LST.filled <- LST.in
#   
#   ptm <- proc.time()
# 
#   for (i in 1:nrow(LST.test))
#   {
#     x <- unlist(LST.test[i,])
#     y <- unlist(CLST[i,])
#     
#     
#     for (l in 1:46)
#     {
#       if(is.na(x[l]))
#         {
#           x[l] <- y[l] 
#         }
#     
#     LST.filled[i,] <- x  
#     }
#   }
# 
#   proc.time() - ptm

# user  system elapsed 
# 1675.05    0.01 1676.25

#__This was explicit gap-filling only for this case
#___leaving it in for the record________________________________________________________________


#   for(i in 1:nrow(LST.in))
#     {
#       if (is.na(LST.in[i,7]))
#       {
#         LST.in[i,7] <-CLST.out[i,7]
#       }
#     }
# 
# 
#   for(i in 1:nrow(LST.in))
#     {
#       if (is.na(LST.in[i,3]))
#       {
#         LST.in[i,3] <-CLST.out[i,3]
#       }
#     }
# 
#   for(i in 1:nrow(LST.in))
#     {
#       if (is.na(LST.in[i,9]))
#       {
#         LST.in[i,9] <-CLST.out[i,9]
#       }
#     }
# 
#   LST.in[99,7] <- 13445.5 #manual fixes for too many NA columns in a row
#   LST.in[100,7] <- 13446.5
#   LST.in[173,7] <- 13476.5
# 
#   LST.in[704,] <- LST.in[703, ] #manual fixes for NA rows
#   LST.in[656,] <- LST.in[655, ]
#   LST.in[608,] <- LST.in[607, ]
# 
#   
# 
################
  
  tLST.in <- t(LST.in)
  tLST.out <- na.spline(tLST.in)
# 
#   na_count <-sapply(LST.in, function(y) sum(length(which(is.na(y)))))
#   na_count <- data.frame(na_count)

  
  CLST <- tLST.out*0.02-273.15
  tCLST <- t(CLST)
  CLST.out <- as.data.frame(tCLST)

########################
#   myFunc <- function(x)
#               {
#               na.spline(x)
#               
#               x
#               }
#               
#   LST.filled <- apply(LST.in, 2, function(x) myFunc(x))

#__These functions may no longer be necessary when using na.spline()_____

#   myFunc2 <- function(x)
#     {
#       if (is.na(x[45]))
#       {
#         x[45] <- x[44]
#       }
#       x
#     }
# 
#   tLST.filled <- t(LST.filled)
#   LST.filled <- apply(tLST.filled, 1, function(x) myFunc2(x))
#   
#   myFunc3 <- function(x)
#     {
#       if (is.na(x[46]))
#       {
#         x[46] <- x[45]
#       }
#       x
#     }
#   
#   tLST.filled <- t(LST.filled)
#   LST.filled <- apply(tLST.filled, 1, function(x) myFunc3(x))
#   tLST.filled <- t(LST.filled)
#   LST.out <- as.data.frame(tLST.filled)
###############################

  CLST.out$PointID <- PointID
  colnames(CLST.out)[1:46] <- colnames(LST.in)[1:46]

  write.dbf(CLST.out, file = paste0("LST", year, "_", medBasin, "_interp.dbf"))

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

  setwd(paste0(mainPath, longBasin))

  out_Preds <- read.dbf(paste0("LST", year, "_", medBasin, "_interp.dbf")) 
  
  setwd(paste0(mainPath, longBasin))
  
  weights <- read.dbf("CW_area_wgts.dbf")
  weights[weights$area_wgts<0.01, "area_wgt"] = 0.01
  
  rcas <- unique(unlist(weights$rca_id))
  rca_zonal <- matrix(ncol=46, nrow=length(rcas))


  l <- 1
  for(i in rcas)	
      {
      pixels <- weights[weights$rca_id == i, "ID"]
      wgts <- weights[weights$rca_id == i, "area_wgt"]
      
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


write.dbf(rca_zonal, file = paste0("LST", year, "_", medBasin, "_RCA.dbf"))

#################################################################
# Logger clean-up and parse
# Filters the larger CBS data set down to only sites within the 
# basins of interest
#################################################################

# setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Clearwater/")
# 
#   Log.in <- read.csv("Clearwater_CBS_logger_data_cleaned.csv")
#   ID.in <- read.csv("Cwater_site_ID.csv")
# 
#   ID.in$Orig_ID <- as.character(ID.in$Orig_ID)
#   SiteID <- unique(ID.in$Orig_ID)
#   SiteID <- as.matrix(SiteID)
#   Log.out <- data.frame (mup = NULL)
# 
#   for (i in SiteID) 
#     { 
#       Log.site <- Log.in[Log.in$OrigID == i,]
#       
#       Log.out <- rbind(Log.out, Log.site)
#       
#     }
# 
#   write.csv(x = Log.out, file = "Clearwater_only_CBS_logger_data_cleaned.csv", row.names = F)

#################################################################
#Logger summary
#Fills in date gaps with NAs and generates 8-day means 
#################################################################

  setwd(paste0(mainPath, longBasin))

  ID.in <- read.dbf(paste0(basin, "_sites_elev.dbf"), as.is=TRUE)
  LST.in <- read.dbf(paste0("LST", year, "_", medBasin, "_RCA.dbf"), as.is=TRUE)

  Log.in <- read.csv("CW_CBS_2012_logger_data.csv", stringsAsFactors = FALSE)
  Log.in$OrigID <- as.character(Log.in$OrigID)


  
  
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    #specifies those dates to be numeric in order to drop the leading zeros
  #newnamesnum <- as.numeric(colnames(LST.in)[1:46]);
  colnames(LST.in)[1:46] <- newnamesnum;
  
  dates <- matrix(nrow=365, ncol=1)
  dates[,1] <- 1:365
  colnames(dates) <- c("Days")
  
  SiteID <- unique(Log.in$OrigID)
  SiteID <- as.matrix(SiteID)
  Log.8Day.out <- data.frame (mup = NULL)
  
    for (i in SiteID) 
      { 
        Log.site <- Log.in[Log.in$OrigID == i,]
        full.year <- merge(dates, Log.site, by.x = "Days", by.y = "JulDay", all.x = TRUE)
        eightday <- rollapply(full.year$DailyMaximum, 8, mean, fill=NA, align = "left")
        eightday <- as.matrix(eightday)
        full.year$Mx8D <- eightday
        full.year$Mx8D[361] <- mean(full.year$DailyMaximum[361:365])
        Log.8Day.out <- rbind(Log.8Day.out, full.year)
        
      }

  write.csv(x = Log.8Day.out, file = paste0(basin, "_", yearPath, "_8Day_logger_data.csv", row.names = F))


#################################################################
#Logger prediction modeling part
#Generates models for both LST-only and LST + Julian Day datasets
#################################################################
 
  Log.8Day.out <- read.csv(paste0(basin, "_", yearPath, "_8Day_logger_data.csv"), stringsAsFactors=FALSE)
  
  setwd(paste0(mainPath, yearPath))
 
  LST.Log.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.8Day.out[Log.8Day.out$OrigID  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$rca_id[ID.in$Orig_ID == i]
      if (length(RCAID) > 0)
      {
      Elev <- ID.in$ELEV[ID.in$Orig_ID == i]
      
      LST.site <- matrix(ncol=3, nrow=46)
      LST.site[,1] <- as.numeric(unlist(colnames(LST.in)[1:46]))
      LST.site[,2] <- unlist(LST.in[RCAID,1:46])
      LST.site <- data.frame(LST.site)
      colnames(LST.site) <- c("JulDay", "LST", "Elev")
      LST.site[3] <- Elev
      LST.Log.site <- merge(LST.site, Log.site, by.x = "JulDay", by.y = "Days", all.x=TRUE, all.y = FALSE)
      
      LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
      }
    }
    
  ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
  NoNA.xyz <- LST.Log.out[ind,]
  
  NoNA.xyz <- NoNA.xyz[,c(13, 2, 1, 3, 6)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  plot(NoNA.xyz$x, NoNA.xyz$y)
  

  write.csv(x=NoNA.xyz, file= paste0(basin, "_", yearPath, "_8Day_model_data.csv"), row.names = FALSE)
            
####################################
# Parsing out the subwatersheds
######################################

  NoNA.xyz <- read.csv(paste0(basin, "_", yearPath, "_8Day_model_data.csv"), stringsAsFactors=FALSE)
  
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:maxrow,]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]

  
  basin <- "Lochsa"
  setwd(paste0(mainPath, longBasin))

  Sub_ID.in <- read.dbf(paste0(basin,"_sites_elev.dbf"), as.is=TRUE)
  
  Sub_IDs <- unique(Sub_ID.in$Orig_ID)

  #Sub_data <- summer[summer$SiteName %in% Sub_IDs, ]
  #Sub_data$y <- log(Sub_data$y)



  Sub_data <- NoNA.xyz[NoNA.xyz$SiteName %in% Sub_IDs, ]
  Sub_data$y <- log(Sub_data$y)
  y <- Sub_data$y
  x <- Sub_data$x
  z <- Sub_data$z
  e <- Sub_data$e
  plot(x, y)
  
  mod <- lm(y ~ x + I(x^2) + e)
  sum_mod <- summary(mod)

  maxrow <- which.max(Sub_data$y)
  data.sp <- Sub_data[1:maxrow,]
  data.fall <- Sub_data[maxrow:nrow(Sub_data),]

  
  
  
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
pred.out[,5] <- "2012"
colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Ent_2012_summer_only.csv",sep = ",", col.names=F)  

plot(pred.out[,1], pred.out[,2])
abline(0,1)

################################
# full year
###################################

  setwd(paste0(mainPath, longBasin, "/", yearPath))

  y <- NoNA.xyz$y
  x <- NoNA.xyz$x
  z <- NoNA.xyz$z
  e <- NoNA.xyz$e
  plot(x, y)
  


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
  pred.out[,5] <- "2012"
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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall_log.csv"), sep = ",", col.names=F)  

  sp_pred <- subset(pred.out, z > 181 & z < 258)
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

  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall_log.csv"),sep = ",", col.names=F)  

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
  
  setwd(paste0(mainPath, longBasin))
  elev.in <- read.csv(paste0(basin, "_rca_elev.csv"))
  LST.in <- read.dbf(paste0("LST", year, "_", medBasin, "_RCA.dbf"), as.is=TRUE)
  



  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));
  colnames(LST.in)[1:46] <- newnamesnum;

  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")

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

  for (i in 1:length(rcas))  
    {
      x <- unlist(LST.sum[i,])
      
      j <- 185
      for (l in 1:12)
      {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[13] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
       j <- j + 8}
      
      LogPred.out[i,1:12] <- x [1:12] 
    }

  LogPred.out <- as.data.frame(LogPred.out)

  #LogPred.out[,1:12] <- exp(LogPred.out[,1:12])

  LogPred.out[LogPred.out< -0.5] = -0.5
  LogPred.out$Basin_RCA <- paste0(basin, "_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:12]))
  varName <- paste0("Tmx_", year)
  names.out <- sprintf("%s_%03d", varName, namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]


write.dbf(LogPred.out, file = paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf")) 


##############
# went back and calc-ed the model metrics for sping & fall together:
###############

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/2012/")

Log.in <- read.csv("Ent_2012_8Day_model_data.csv")
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

############################
# 15Jul-31Aug ------------------------------------------------------------------
#############################

  setwd(paste0(mainPath, longBasin))
  basin <- "CW"
  popn <- "STHD"
  varName <- paste0(basin, "_", popn)

  Au.in <- read.csv(paste0(basin, "_RCA_AUs_IP.csv"), stringsAsFactors = FALSE)
  
 
  basin <- "Lochsa"

  ST_rcas <- Au.in[Au.in$basin == "Steelhead",]
  setwd(paste0(mainPath, longBasin, "/", yearPath))


  
  Max.in <- read.dbf(paste0("predt", yearPath, "_", basin, "_8D_Max_summer.dbf"))
  colnames(Max.in)<-gsub("Tmx_12_", "", colnames(Max.in))

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

#   Max.rep <- Max.in[,c(3,16,4,16,5,16,6,16,7,16,8,16,9)]
#   toFill <- t(Max.rep)
#   filled <- t(interpNA(toFill, method="linear"))
#   dim(filled)
# 
#   toFill <- cbind(filled, Max.in[,16])
#   toFill.rep <- toFill[,c(1,14,2,14,3,14,4,14,5,14,6,14,7,14,8,14,9,14,10,14,11,14,12,14,13)]
#   ttoFill <- t(toFill.rep)
#   filled <- t(interpNA(ttoFill, method="linear"))
#   dim(filled)
# 
#   toFill <- cbind(filled, Max.in[,16])
#   toFill.rep <- toFill[,c(1,26,2,26,3,26,4,26,5,26,6,26,7,26,8,26,9,26,10,26,11,26,12,26,13,26,14,26,15,26,16,26,17,26,18,26,19,26,20,26,21,26,22,26,23,26,24,26,25)]
#   ttoFill <- t(toFill.rep)
#   filled <- t(interpNA(ttoFill, method="linear"))
#   dim(filled)
# 
#   class(filled) <- "numeric"
#   filled.out <- as.data.frame(filled, stringsAsFactors=FALSE)
# 
#   filled.out$RCAID <- Max.in$RCAID
#   colnames(filled.out)[1:49] <- 201:249

  Max.au <- merge(Max.out, Au.in, by.x = "RCAID", by.y = "rca_id")




#_________all_______________________  
rcas <- unique(Max.au$RCAID)

SumSumm.out <- data.frame ("RCAID" = rcas)


# count of days in exceedence -----------------------------------------


  for (i in 1:length(rcas)) 
    { 
      l <- rcas[i]
      MaxRCA <- Max.au[Max.au$RCAID == l,] #grab days for one RCA 
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


  Aus <- unique(Max.au$Lochsa_AU)
  Aus <- Aus[-1]
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
      SumAU <- SumSummAU[SumSummAU$Lochsa_AU == Aus[i] & SumSummAU$Lochsa == "Steelhead",]
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


########
#________chinook_______________________
########

Aus <- unique(Max.au$AU_ID)
rcas <- unique(CH_rcas$rca_id)

########
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
write.csv(SumAU.out, file = "Ent_AU_2012_21Jul_31Aug_ExPcnt_summary_CHNK_IP.csv", row.names = F)


#######
#_________________pie charts_____________
#######
`

library(plotrix)

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/"
yearPath <- "2012"
setwd(paste0(mainPath, longBasin, "/", yearPath))
graphPath <- "D:/OneDrive/work/research/CHaMP/graphics/EP/"
basin <- "Lochsa"
popn <- "STHD"


  SumAU.out <- read.csv(paste0(basin, "AU_", yearPath, "_21Jul_31Aug_ExPcnt_summary_", popn, "_IP.csv"), stringsAsFactors = FALSE)
  SumAU.out$AU_Code <- as.character(SumAU.out$AU_Code)
  Aus <- SumAU.out$AU_Code


  metric <- "20"
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