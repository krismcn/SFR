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
mainPath <- "D:/OneDrive/work/research/CHaMP/"
basin <- "BNG"
longBasin <- "Big_Navarro"
subDir <- "GIS/LST/LST_s3"
setwd(paste0(mainPath, "/", subDir))

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

  LST.in <- read.dbf("BNG_13_LST.dbf")
  GrPolID <- LST.in[,1]
  LST.in <- LST.in[,2:47]
  LST.in[LST.in<0.1] = NA
  
  LST.names<-colnames(LST.in)
  LST.names <- order(LST.names)
  LST.in<- LST.in[,(LST.names)]
  
  tLST.in <- t(LST.in)
  tLST.out <- na.spline(tLST.in)
  CLST <- tLST.out*0.02-273.15
  tCLST <- t(CLST)
  CLST.out <- as.data.frame(tCLST)
  
  CLST.out$GRID_CODE <- GrPolID
  
  colnames(CLST.out)[1:46] <- LST.names

  subDir <- "GIS/coverages/Big_Navarro_Garcia"
  setwd(paste0(mainPath, "/", subDir))  

  write.dbf(CLST.out, file = "LST13_BNG_interp.dbf")

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

out_Preds <- read.dbf("LST13_BNG_interp.dbf") #use the appropriate read statement

weights <- read.dbf("BNG_rca_area_wgts.dbf")
weights[weights$area_wgts<0.01, "area_wgts"] = 0.01
weights[weights$area_wgts<0.01, "area_wgt"] = 0.01

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

write.dbf(rca_zonal, file = "LST13_BNG_RCA.dbf")



#################################################################
#Logger prediction modeling part
# First section reads in the Excel file and processes the sheets
#################################################################

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  basin <- "BNG"
  longBasin <- "Big_Navarro"
  subDir <- "CA_Temperature_data"
  setwd(paste0(mainPath, longBasin, "/", subDir))
  fileName <- "Complete_Hobo_Data_2012-Present.xlsx"
  #fileName <- "Mendocino_temperature_data2.xlsx"

  sheets <- excel_sheets(paste0(fileName))
  merge.out <- data.frame(SiteID=NULL, yrJul=NULL, Mean=NULL, Max=NULL, obsCount=NULL, countFlag=NULL, meanFlag=NULL)
  out.stats <- data.frame(SiteID=Null, FirstDay=NULL, LastDay=NULL)

  for (i in 2:length(sheets))
  {
    #       coltypes <- c("numeric", "date", "date", "text", "numeric")
    #       colnames <- c("Number", "Date", "Time", "12H", "Temp")
    #       Site.data <- read_excel(paste0(fileName), sheet = i, col_names =  colnames, col_types = coltypes, na = "", skip = 2)
    #       ind <- apply(Site.data, 1, function(x) !any(is.na(x)))
    #       good.data <- Site.data[ind,]
    #       good.data <- good.data[good.data$QA == "pass",]
    #       
    ####### if Excel is being weird#######
    Site.data <- read_excel(paste0(fileName), sheet = i, col_names = FALSE, skip = 2)
    colnames(Site.data) <- c("Number", "Date", "Time", "12 hour", "Temp")
    Site.data$Date <- as.POSIXlt(Site.data$Date)
    Site.data$Time <- as.POSIXct(Site.data$Time)
    Site.data$Temp <- as.numeric(Site.data$Temp)
    good.data <- Site.data
    ########################################
    
    good.data$JulDay <- strptime(good.data$Date, "%F")$yday+1
    good.data$year <- as.numeric(format(good.data$Date, '%y'))
    good.data$yrJul <- as.numeric(paste(good.data$year, good.data$JulDay, sep=""))
    
    firstDay <- good.data$yrJul[1]
    num_obs <- as.integer(dim(good.data)[1])
    lastDay <- good.data$yrJul[num_obs]
    num_days <- length(unique(good.data$yrJul))
    unique.Days <- unique(good.data$yrJul)
    Site.stats <- data.frame(sheets[i], firstDay, lastDay)
    out.stats <- rbind(out.stats, Site.stats)
    
    out.summary <- data.frame(SiteID=numeric(num_days), yrJul=numeric(num_days), Mean=numeric(num_days), Max=numeric(num_days), obsCount=numeric(num_days), countFlag=numeric(num_days), meanFlag=numeric(num_days))
    out.summary$SiteID <- sheets[i]
    out.summary$yrJul <- unique(good.data$yrJul)
    
    
    for (j in unique.Days)
      {
        day_indeces <- which(good.data$yrJul == j)
        temps_at_day <- good.data[day_indeces,]
        out.summary$obsCount[out.summary$yrJul == j] <- length(day_indeces)
        out.summary$countFlag[out.summary$yrJul == j] <- as.character(length(day_indeces) == 24) #flags days with <> 24 datapoints
        out.summary$Mean[out.summary$yrJul == j] <- mean(temps_at_day$Temp)
        out.summary$Max[out.summary$yrJul == j] <- max(temps_at_day$Temp)
        out.summary$meanFlag[out.summary$yrJul == j] <- head(as.character(out.summary$Mean[out.summary$yrJul == j] > 32),1) #flags for days with a mean temp >32'
        filename.data <-paste0("site_data_",as.character(sheets[i]),".csv")
        write.table (x=out.summary,append=F,row.names=F,file=filename.data,sep = ",", col.names=T)  #writes out each sites data as an individual file
      }
    
    merge.out <- rbind(merge.out, out.summary)
  }


write.table(x=merge.out, append = F, file= paste0(basin, "_logger_data.csv"),sep = ",", col.names=T, row.names = FALSE)
write.table(x=out.stats, append = F, file= paste0(basin, "_site_date_ranges.csv"),sep = ",", col.names=T, row.names = FALSE)

#################################################################
#Logger summary
#Fills in date gaps with NAs and generates 8-day means 
#################################################################

  setwd(paste0(mainPath, longBasin))
  
  
  dates <- matrix(nrow=365, ncol=1)
  dates[,1] <- 1:365
  colnames(dates) <- c("Days")
  
  SiteID <- unique(merge.out$SiteID)
  SiteID <- as.matrix(SiteID)
  Log.data<- merge.out[,1:4]
  Log.data$JulDay <- as.numeric(substring(Log.data$yrJul,3,5))
  Log.data$year <- as.numeric(substring(Log.data$yrJul,1,2))
  year <- unique(Log.data$year)
  year <- year[-5]
  

  for (j in year)
    {
      Log.in <- Log.data[Log.data$year == j,]
      Log.8Day.out <- data.frame (mup = NULL)
      yearPath <- j
      
      for (i in SiteID) 
        { 
          Log.site <- Log.in[Log.in$SiteID == i,]
          if (dim(Log.site)[1] > 8)
            {
              Log.site <- Log.in[Log.in$SiteID == i,]
              full.year <- merge(dates, Log.site, by.x = "Days", by.y = "JulDay", all.x = TRUE)
              eightday <- rollapply(full.year$Mean, 8, mean, fill=NA, align = "left")
              eightday <- as.matrix(eightday)
              full.year$Mx8D <- eightday
              full.year$Mx8D[361] <- mean(full.year$Mean[361:365])
              Log.8Day.out <- rbind(Log.8Day.out, full.year)
              #               windows()
              #               plot(full.year$Days, full.year$Mean, main = c(i,j))
              remove(Log.site)
            } else {
              remove(Log.site)
            }
        }
      write.csv(x = Log.8Day.out, file = paste0(basin, "_", yearPath, "_Mn_8Day_logger_data.csv"), row.names = F)
    }

##########################################################

  setwd(paste0(mainPath, longBasin))
  
  yearPath <- 13
  Log.8Day.out <- read.csv(file = paste0(basin, "_", yearPath, "_Mn_8Day_logger_data.csv"))
  Logger <- na.omit(Log.8Day.out)
  Logger$yrJul <- as.numeric(sprintf("%1d%03d", Logger$year, Logger$Days))
  write.csv(x = Logger, file = paste0(basin, "_", yearPath, "_all_sites_Mn_8Day_logger_data.csv"), row.names = F)
  
  mainPath <- "D:/OneDrive/work/research/CHaMP/"
  basin <- "BNG"
  longBasin <- "Big_Navarro"
  subDir <- "CHaMP_data"
  setwd(paste0(mainPath, subDir, "/", longBasin))

  Log.in <- read.csv(file = paste0(basin, "_", yearPath, "_all_sites_Mn_8Day_logger_data.csv"))
  ID.in <- read.csv("BNG_sites_elev.csv")

  longBasin <- "Big_Navarro_Garcia"
  subDir <- "GIS/coverages"
  setwd(paste0(mainPath, subDir, "/", longBasin))
  LST.in <- read.dbf("LST13_BNG_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  
  
  longBasin <- "Big_Navarro"
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  yearPath <- "2013"
  setwd(paste0(mainPath, longBasin, "/", yearPath))

  
  Log.in$SiteID <- as.character(Log.in$SiteID)
  Log.in$SiteID<-gsub("GRTS#", "", Log.in$SiteID)
  SiteID <- unique(Log.in$SiteID)
  SiteID <- as.matrix(SiteID)
  ID.in$SiteID <- as.character(ID.in$GRTS)
  colnames(Log.in)[1] <- "JulDay"
  LST.Log.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$SiteID  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$rca_id[ID.in$SiteID == i]
      Elev <- ID.in$DEM_10m[ID.in$SiteID == i]
      
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

  NoNA.xyz <- NoNA.xyz[,c(6, 2, 1, 3, 4)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  plot(NoNA.xyz$x, NoNA.xyz$y)

  write.csv(x=NoNA.xyz, file=paste0(basin, "_", yearPath, "_8Day_model_data.csv"), row.names = FALSE)

  NoNA.xyz <- orderBy(~z, NoNA.xyz)

  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:maxrow,]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]



################################
# full year
###################################

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
  pred.out[,5] <- yearPath
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_full_year.csv"),sep = ",", col.names=F)  
  
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
  plot(x, y)

  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y<0] = 0.0
  plot(y, pred.y)


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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F)  

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
  points(y, pred.y, col="blue", pch=16)
  abline(0,1)


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

  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall.csv"),sep = ",", col.names=F) 

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

write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", yearPath, "_mod_coeffs_Mn.csv"), sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", yearPath, "_mod_metrics_Mn.csv"), sep = ",", col.names=T)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################

  mainPath <- "D:/OneDrive/work/research/CHaMP/"
  yearPath <- "2013"  
  longBasin <- "Big_Navarro_Garcia"
  subDir <- "GIS/coverages"
  setwd(paste0(mainPath, subDir, "/", longBasin))
  LST.in <- read.dbf("LST13_BNG_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  elev.in <- read.csv("BNG_RCA_Elev.csv")

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  longBasin <- "Big_Navarro"
  
  setwd(paste0(mainPath, longBasin, "/", yearPath))

  colnames(LST.in)[1:46] <- newnamesnum
  
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")
  
 
  coeffs.in <- read.csv(paste0("All_data_", yearPath, "_mod_coeffs_Mn.csv"))
  LogPred.out <- LST.elev
  LogPred.out[,2:48] <- 0
  rcas <- unique(elev.in$rca_id)
  
  
  for (i in 1:length(rcas))  
      {
        x <- unlist(LST.elev[i,])
        maxrow <- as.numeric(which.max(x[2:47]))
        midrow <- maxrow + 1
        j <- 1
        day <- as.numeric(colnames(LST.elev)[maxrow])
        
        
        for (l in 2:maxrow)
          {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[48] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
           j <- j + 8}
        k <- midrow
        for (l in midrow:47)     
          {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[48] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
           k <- k + 8}
        fill <- rollmean(x[2:47], 5, align = "left")
        x[(midrow-4):(midrow+4)] <- fill[(maxrow-4):(maxrow+4)]
      LogPred.out[i,2:47] <- x [2:47] 
    }

  LogPred.out <- as.data.frame(LogPred.out)

LogPred.out[LogPred.out< -0.5] = -0.5
LogPred.out$Basin_RCA <- paste0("BNG_", LogPred.out$RCAID)
namesnum <- as.numeric(colnames(LogPred.out[1:12]))

names.out <- sprintf("Tmx_13_%03d", namesnum)
colnames(LogPred.out)[1:12] <- names.out[1:12]


write.dbf(LogPred.out, file = "predt2013_Ent_8D_Max_summer.dbf") 


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

#_________________pie charts_____________


library(plotrix)

###########################################
# Set parameters for each year, pop'n, and metric
###########################################

  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/Entiat/"
  yearPath <- "2013"
  setwd(paste0(mainPath, yearPath))
  
  metric <- "22"
  popn <- "CHNK"
  varName <- paste0("v", metric, "DMax")

  SumAU.out <- read.csv(paste0("Ent_AU_2013_21Jul_31Aug_ExPcnt_summary_", popn, "_IP.csv"))

############################################

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
