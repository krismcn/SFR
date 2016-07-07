############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# stream temperature for a year using MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# Currently set up for the John Day which has 12039 grid cells and 5532 RCAs (those parameters can be changed for other regions).
# Should do a search-and-replace for the output folder (usually dated), and the year being processed, in both a YYYY and a _YY format.

# Edited Aug 2015 to add the PRESS stastic output
# Should look at PRESS output to screen to chose model - only that model will output resids and preds.


          
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
library(readxl)
library(lubridate)
library(zoo)

mainPath <- "D:/OneDrive/work/research/CHaMP/"
basin <- "BNG"
longBasin <- "Big_Navarro"
subDir <- "GIS/LST/LST_s3_2015"
setwd(paste0(mainPath, "/", subDir))

  LST.in <- read.dbf("BNG_15_LST.dbf")
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

  subDir <- "GIS/coverages/Big_Navarro_Garcia"
  setwd(paste0(mainPath, "/", subDir))  

  write.dbf(CLST.out, file = "LST15_BNG_interp.dbf")

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

  out_Preds <- read.dbf("LST15_BNG_interp.dbf") #use the appropriate read statement
  
  weights <- read.dbf("BNG_rca_area_wgts.dbf")
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

  write.dbf(rca_zonal, file = "LST15_BNG_RCA.dbf")
  


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
              full.year$Mn8D <- eightday
              full.year$Mn8D[361] <- mean(full.year$Mean[361:365])
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

  mainPath <- "D:/OneDrive/work/research/CHaMP/"
  longBasin <- "Big_Navarro_Garcia"
  subDir <- "GIS/coverages"
  setwd(paste0(mainPath, subDir, "/", longBasin))

  LST.in <- read.dbf("LST15_BNG_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  ID.in <- read.csv("BNG_sites_elev.csv")  
  HUC <- read.csv("BNG_RCA_HUC10.csv", stringsAsFactors = FALSE)

  
  basin <- "BNG"
  longBasin <- "Big_Navarro"
  subDir <- "CHaMP_data"
  yearPath <- 2015
  yrPath <- 15
  setwd(paste0(mainPath, subDir, "/", longBasin, "/", yearPath))

########################################################

  Log.8Day.out <- read.csv(file = paste0(basin, "_", yrPath, "_Mn_8Day_logger_data.csv"))
  Logger <- na.omit(Log.8Day.out)
  Logger$yrJul <- as.numeric(sprintf("%1d%03d", Logger$year, Logger$Days))
  write.csv(x = Logger, file = paste0(basin, "_", yrPath, "_Mn_8Day_logger_data.csv"), row.names = F)

###########################################################

  Log.in <- read.csv(file = paste0(basin, "_", yrPath, "_Mn_8Day_logger_data.csv"), stringsAsFactors = FALSE)
  
  colnames(Log.in)[1] <- "JulDay"
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  
  
  Log.in$SiteID<-gsub("GRTS#", "", Log.in$SiteID)
  Log.in$SiteID<-gsub("GRST#", "", Log.in$SiteID)
  SiteID <- unique(Log.in$SiteID)
  SiteID <- as.matrix(SiteID)
  ID.in$SiteID <- as.character(ID.in$GRTS_)
  ID_HUC <- merge(ID.in, HUC, by.x = "rca_id", by.y = "rca_id")
  ID_HUC$HUC_10 <- as.character(ID_HUC$HUC_10)
  Log.HUC <- merge(Log.in, ID_HUC, by.x= "SiteID", by.y = "SiteID")

  LST.Log.out <- data.frame (mup = NULL)


  for (i in SiteID) 
    { 
      Log.site <- Log.HUC[Log.HUC$SiteID  == i,]
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

  NoNA.xyz <- NoNA.xyz[,c(12, 2, 1, 3, 4, 21)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName", "HUC")
  plot(NoNA.xyz$x, NoNA.xyz$y)
  
  HUCs <- unique(NoNA.xyz$HUC)

  write.csv(x=NoNA.xyz, file=paste0(basin, "_", yearPath, "_8Day_Mn_model_data.csv"), row.names = FALSE)

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
  plot(z, y)

  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y<0] = 0.0
  plot(pred.y, y)

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
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_Mn_full_year.csv"),sep = ",", col.names=F)  
  
  plot(pred.out[,1], pred.out[,2])
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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_Mn_sp_fall.csv"),sep = ",", col.names=F)  

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
  pred.out[,5] <- yearPath

  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_Mn_sp_fall.csv"),sep = ",", col.names=F) 


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


######################################
# by HUC
####################################

#   m <- HUCs[5]
#   data.huc <- NoNA.xyz[NoNA.xyz$HUC == m,]
#   length(unique(data.huc$SiteName))
#   dim(data.huc)
#   m


  HUCs <- unique(NoNA.xyz$HUC)
  coeffs.huc <- data.frame (mup = NULL)
  metrics.huc <- data.frame (mup = NULL)
  preds.huc <- data.frame (mup = NULL)

    for(m in HUCs)
      {
        
        data.huc <- NoNA.xyz[NoNA.xyz$HUC == m,]
        data.huc <- orderBy(~z, data.huc)
        maxrow <- which.max(data.huc$y)
        midrow <- maxrow-1
        data.sp <- data.huc[1:midrow,]
        data.fall <- data.huc[maxrow:nrow(data.huc),]
        
        
        coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bJul=numeric(2), bElev=numeric(2))
        metrics_out <- data.frame(PRESS=numeric(2), p2=numeric(2), r2=numeric(2), RMSEP=numeric(2), RMSE=numeric(2))
        rownames(metrics_out) <- c("Spring", m)
        rownames(coeffs_out) <- c("Spring", m)
        
        
        y <- data.sp$y
        x <- data.sp$x
        z <- data.sp$z
        e <- data.sp$e
        plot(x, y, main=m)
        plot(z, y, main=m)
        mod <- lm(y ~ x + I(x^2) + z + e)
        sum_mod <- summary(mod)
        coeffs <- as.matrix(coefficients(mod))
        pred.y <- predict(mod)
        pred.y[pred.y<0] = 0.0
        plot(y, pred.y, main=m)
        
        
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
        #colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
        write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_", m, "_Mn_sp_fall.csv"),sep = ",", col.names=F)  
        preds.huc <- rbind(preds.huc, pred.out)
        
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
        plot(x,y, main=m)
        plot(z,y, main=m)
        
        mod <- lm(y ~ x + I(x^2) + z + e)
        sum_mod <- summary(mod)
        coeffs <- as.matrix(coefficients(mod))
        pred.y <- predict(mod)
        pred.y[pred.y<0] = 0.0
        plot(y, pred.y, col="blue", pch=16, main=m)
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
        pred.out[,5] <- yearPath
        preds.huc <- rbind(preds.huc, pred.out)
        
        write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_", m, "_Mx_sp_fall.csv"),sep = ",", col.names=F) 
        
        
        
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
        
        write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", yearPath, "_", m, "_mod_coeffs_Mn.csv"), sep = ",", col.names=T)
        
        write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", yearPath, "_", m, "_mod_metrics_Mn.csv"), sep = ",", col.names=T)
        
        coeffs.huc <- rbind(coeffs.huc, coeffs_out)
        metrics.huc <- rbind(metrics.huc, metrics_out)
      }

write.csv(x=coeffs.huc, row.names=T, file = paste0("All_HUCs_", yrPath, "_mod_coeffs_Mn.csv"))

write.csv(x=metrics.huc, row.names=T, file = paste0("All_HUCs_", yrPath, "_mod_metrics_Mn.csv"))
########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day Max temp estimates by HUC
########################################################################################################
  
  mainPath <- "D:/OneDrive/work/research/CHaMP/"
  yearPath <- "2015"  
  longBasin <- "Big_Navarro_Garcia"
  subDir <- "GIS/coverages"
  setwd(paste0(mainPath, subDir, "/", longBasin))
  LST.in <- read.dbf("LST15_BNG_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  elev.in <- read.csv("BNG_RCA_Elev.csv")
  HUC.in <- read.csv("BNG_RCA_HUC10.csv", stringsAsFactors = FALSE)

  basin <- "BNG"
  longBasin <- "Big_Navarro"
  subDir <- "CHaMP_data"
  yearPath <- 2015
  yrPath <- 15
  setwd(paste0(mainPath, subDir, "/", longBasin, "/", yearPath))
  
  coeffs.in <- read.csv(paste0("All_HUCs_", yrPath, "_mod_coeffs_Mn.csv"))
  colnames(LST.in)[1:46] <- newnamesnum
  
  elev_HUC <- merge(elev.in, HUC.in, by.x = "rca_id", by.y = "rca_id")
  LST.elev <- merge(LST.in, elev_HUC, by.x = "RCAID", by.y = "rca_id")
  
  
  rcas <- unique(elev.in$rca_id)
  
  
  LogPred.out <- LST.elev
  LogPred.out[,2:47] <- 0
  
    for(i in 1:length(rcas))
      {
        x <- LST.elev[i,2:49]
        
        maxrow <-as.numeric(which.max(x[1:46]))
        day <- as.numeric(colnames(LST.elev)[maxrow])
        midrow <- maxrow - 1
        j <- 1
        
        if (x[48] == 1501010806)
        {
          for (l in 1:midrow)
          {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[47] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
           j <- j + 8}
          k <- day
          for (l in maxrow:46)     
          {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[47] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
           k <- k + 8}
          fill <- rollmean(unlist(x[1:46]), 3, align = "left")
          x[(midrow-4):(midrow+4)] <- fill[(midrow-4):(midrow+4)]
        } else
          if (x[48] == "ocean")
          {
            for (l in 1:midrow)
            {x[l] <- x[l] * coeffs.in$bLST[3] + j * coeffs.in$bJul[3] +  x[47] * coeffs.in$bElev[3] +  x[l]^2 * coeffs.in$bLST2[3] + coeffs.in$Int[3]
             j <- j + 8}
            k <- day
            for (l in maxrow:46)     
            {x[l] <- x[l] * coeffs.in$bLST[4] + k * coeffs.in$bJul[4] +  x[47] * coeffs.in$bElev[4] + x[l]^2 * coeffs.in$bLST2[4] + coeffs.in$Int[4]
             k <- k + 8}
            fill <- rollmean(unlist(x[1:46]), 3, align = "left")
            x[(midrow-4):(midrow+4)] <- fill[(midrow-4):(midrow+4)]
          } else
            if (x[48] == 1801010809)
            {
              for (l in 1:midrow)
              {x[l] <- x[l] * coeffs.in$bLST[5] + j * coeffs.in$bJul[5] +  x[47] * coeffs.in$bElev[5] + x[l]^2 * coeffs.in$bLST2[5] + coeffs.in$Int[5]
               j <- j + 8}
              k <- day
              for (l in maxrow:46)     
              {x[l] <- x[l] * coeffs.in$bLST[6] + k * coeffs.in$bJul[6] +  x[47] * coeffs.in$bElev[6] + x[l]^2 * coeffs.in$bLST2[6] + coeffs.in$Int[6]
               k <- k + 8}
              fill <- rollmean(unlist(x[1:46]), 3, align = "left")
              x[(midrow-4):(midrow+4)] <- fill[(midrow-4):(midrow+4)]
            } else
              if (x[48] == 1801010808)
              {
                for (l in 1:midrow)
                {x[l] <- x[l] * coeffs.in$bLST[7] + j * coeffs.in$bJul[7] + x[47] * coeffs.in$bElev[7] + x[l]^2 * coeffs.in$bLST2[7] + coeffs.in$Int[7]
                 j <- j + 8}
                k <- day
                for (l in maxrow:46)     
                {x[l] <- x[l] * coeffs.in$bLST[8] + k * coeffs.in$bJul[8] +  x[47] * coeffs.in$bElev[8] + x[l]^2 * coeffs.in$bLST2[8] + coeffs.in$Int[8]
                 k <- k + 8}
                fill <- rollmean(unlist(x[1:46]), 3, align = "left")
                x[(midrow-4):(midrow+4)] <- fill[(midrow-4):(midrow+4)]
              }else
                if (x[48] == 1801010807)
                {
                  for (l in 1:midrow)
                  {x[l] <- x[l] * coeffs.in$bLST[9] + j * coeffs.in$bJul[9] +  x[47] * coeffs.in$bElev[9] + x[l]^2 * coeffs.in$bLST2[9] + coeffs.in$Int[9]
                   j <- j + 8}
                  k <- day
                  for (l in maxrow:46)     
                  {x[l] <- x[l] * coeffs.in$bLST[10] + k * coeffs.in$bJul[10] +  x[47] * coeffs.in$bElev[10] + x[l]^2 * coeffs.in$bLST2[10] + coeffs.in$Int[10]
                   k <- k + 8}
                  filled <- rollmean(unlist(x[1:46]), 3, align = "left")
                  x[(midrow-4):(midrow+4)] <- filled[(midrow-4):(midrow+4)]
                }else
                  if (x[48] == "")
                  {
                    x[1:46] <- NA
                  }
        LogPred.out [i,2:49] <- x
      }
    
  
  plot(1:46, LogPred.out[LogPred.out$RCAID == 173,2:47])
  points(1:46, LogPred.out[LogPred.out$RCAID == 925,2:47], pch=16, col="blue")
  points(1:46, LogPred.out[LogPred.out$RCAID == 371,2:47], pch=16, col="green")
  points(1:46, LogPred.out[LogPred.out$RCAID == 599,2:47], pch=16, col="red")
  points(1:46, LogPred.out[LogPred.out$RCAID == 1012,2:47], pch=16, col="lightblue")
  points(1:46, LogPred.out[LogPred.out$RCAID == 307,2:47], pch=16, col="purple")
  
  plot(newnamesnum[22:34], LogPred.out[LogPred.out$RCAID == 173,23:35], ylim=c(10,25))
  points(newnamesnum[22:34], LogPred.out[LogPred.out$RCAID == 925,23:35], pch=16, col="blue")
  points(newnamesnum[22:34], LogPred.out[LogPred.out$RCAID == 371,23:35], pch=16, col="green")
  points(newnamesnum[22:34], LogPred.out[LogPred.out$RCAID == 599,23:35], pch=16, col="red")
  points(newnamesnum[22:34], LogPred.out[LogPred.out$RCAID == 1012,23:35], pch=16, col="lightblue")
  points(newnamesnum[22:34], LogPred.out[LogPred.out$RCAID == 307,23:35], pch=16, col="purple")

  LogPred.out[LogPred.out<0] = 0.0
  LogPred.out$Basin_RCA <- paste0("BNG_", LogPred.out$RCAID)
  
  names.out <- sprintf("Tmn_15_%03d", newnamesnum)
  colnames(LogPred.out)[2:47] <- names.out[1:46]
  
  LogPred.out <- LogPred.out[,-48]
  write.dbf(LogPred.out, file = "predt2015_BNG_8D_Mn.dbf") 


#________________________________________________________
# summary for the summer 
# 26June-22Sept ------------------------------------------------------------------



  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  yearPath <- "2015"  
  longBasin <- "Big_Navarro"
  
  library(timeSeries)
  library(lattice)
  library(foreign)
  library(doBy)
  library(qpcR)
  library(pls)
  library(boot)
  library(Hmisc)
  
  setwd(paste0(mainPath, longBasin, "/", yearPath))
#####
#_________Mean_______________________ 
#####
Mean.in <- read.dbf("predt2015_BNG_8D_Mn.dbf") 

rcas <- unique(Mean.in$RCAID)

SumSumm.out <- data.frame ("RCAID" = rcas)

  for (i in 1:length(rcas)) 
    { 
      l <- rcas[i]
      MeanRCA <- Mean.in[Mean.in$RCAID == l,] #grab days for one RCA 
      MaxMean <- max(MeanRCA[2:47])
      SDMn <- sd(MeanRCA[2:47])
      MnMnSummer <- mean(unlist(MeanRCA[24:35]))
      SumSumm.out$MxMn[i] <- MaxMean 
      SumSumm.out$SdMn[i] <- SDMn
      SumSumm.out$MnMnSummer[i] <- MnMnSummer
    } 



colnames(SumSumm.out) <- c("RCAID", "AnMxMn14","AnSDMn14", "SuMnMn14")

write.dbf(SumSumm.out, file = "BNG_2015_26June_22Sept_mean_summary_All.dbf")
write.csv(SumSumm.out, file = "BNG_2015_26June_22Sept_mean_summary_All.csv", row.names = F)
