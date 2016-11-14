############################################################################################################
# This set of R scripts combines a the series of scripts that make up the work flow to predict and validate 
# minimum stream temperature for a year using MODIS 1km LST data
# The initial input is two LST_YY.dbfs which are tables (columns = days, rows = grid cells) with the grid value for each grid cell in the spatial extent
# 
# Edited Aug 2014 to add the PRESS stastic output
# Edited 10 March 2016 to process calendar year data into winter year data and add more model diagnostics
# 


          
##############################################################################################
# This section reads in the 1km LST data for a year, splices the second half of year1 to the first half of year2
# then fills any gaps using a spline function
##############################################################################################

mainPath <- "D:/OneDrive/work/research/CHaMP/"
basin <- "Ent"
longBasin <- "Entiat"
subDir <- "CHaMP_data/"
modelDir <- "Min_models/"
yrPath1 <- 14
yrPath2 <- 15
yearPath1 <- 2014
yearPath2 <- 2015
setwd(paste0(mainPath, "/", subDir, "/", longBasin))

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

  LST.in.1 <- read.csv(paste0("Ent_LST_", yearPath1, ".csv"), stringsAsFactors = FALSE)
  LST.in.2 <- read.csv(paste0("Ent_LST_", yearPath2, ".csv"), stringsAsFactors = FALSE)
  
  colnames(LST.in.1)<-gsub("lsts1", "", colnames(LST.in.1))
  colnames(LST.in.2)<-gsub("X", "", colnames(LST.in.2))
  charnames <- as.numeric(colnames(LST.in.2)[2:47])
  colnames(LST.in.2)[2:47] <- as.numeric(sprintf("%1d%03d", yrPath2, charnames))
  
  LST.in <- cbind(LST.in.1[,27:48], LST.in.2[,3:26]) #takes julday 193-361 from yr1 and appends 1-185 from yr2
                  
  GrPolID <- LST.in.1[,1]
  LST.in[LST.in<0.1] = NA
  
  LST.names<-colnames(LST.in)
  
  tLST.in <- t(LST.in)
  tLST.out <- na.spline(tLST.in)
  CLST <- tLST.out*0.02-273.15
  tCLST <- t(CLST)
  CLST.out <- as.data.frame(tCLST)

  plot(1:46, CLST.out[10,])
  
  CLST.out$GRID_CODE <- GrPolID
  
  colnames(CLST.out)[1:46] <- LST.names

  

  write.dbf(CLST.out, file = paste0("LST_", yrPath1, "_", yrPath2, "_", basin, "_interp.dbf"))

  rm(CLST, CLST.out, LST.in, LST.in.1, LST.in.2, tCLST, tLST.in, tLST.out)

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

  out_Preds <- read.dbf(paste0("LST_", yrPath1, "_", yrPath2, "_", basin, "_interp.dbf")) 
  
  weights <- read.dbf(paste0(basin, "_rca_area_wgts.dbf"))
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
  
  plot(1:46, rca_zonal[1,1:46])
  
  write.dbf(rca_zonal, file = paste0("LST_", yrPath1, "_", yrPath2, "_", basin, "_RCA.dbf"))



  rm(out_Preds, rca_zonal, weights)

##########################################################


  Log.in.1 <- read.csv(file = paste0(basin, "_", yearPath1, "_logger_data.csv"), stringsAsFactors = FALSE)

  Log.in.1$yrJul <- as.numeric(sprintf("%1d%03d", yrPath1, Log.in.1$JulDay))

  Log.in.2 <- read.csv(file = paste0(basin, "_", yearPath2, "_logger_data.csv"), stringsAsFactors = FALSE)

  Log.in.2$yrJul <- as.numeric(sprintf("%1d%03d", yrPath2, Log.in.2$JulDay))
  
  Log.in <- rbind(Log.in.1, Log.in.2)


  LST.in <- read.dbf(paste0("LST_", yrPath1, "_", yrPath2, "_", basin, "_RCA.dbf"))
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))

  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  ID.in <- read.csv(paste0(basin, "_sites_elev.csv"), stringsAsFactors = FALSE)
  
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

  NoNA.xyz <- NoNA.xyz[,c(17, 2, 1, 3, 7)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  plot(NoNA.xyz$x, NoNA.xyz$y, main="LST vs. Logger", xlab="LST", ylab="Logger")

setwd(paste0(mainPath, subDir, longBasin, "/", modelDir, yrPath1, "_", yrPath2, "_Min/"))

  write.csv(x=NoNA.xyz, file=paste0(basin, "_", yrPath1, "_", yrPath2, "_8Day_Min_model_data_LST.csv"), row.names = FALSE)

  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  plot(NoNA.xyz$y)
  plot(NoNA.xyz$x, NoNA.xyz$y)
  
  
#####################################
# cooling/warming legs
######################################
  
  NoNA.xyz <- read.csv(paste0(basin, "_", yrPath1, "_", yrPath2, "_8Day_Min_model_data_LST.csv"), stringsAsFactors=FALSE)
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  
  minrow <- which.min(NoNA.xyz$y)
  split <- NoNA.xyz[minrow, 'z']
  data.cool <- NoNA.xyz[1:minrow,]
  maxrow <- which.max(data.cool$y)
  data.cool <- data.cool[(maxrow+1):dim(data.cool)[1],]
  min <- data.cool$z[1]

  data.warm <- NoNA.xyz[(minrow+1):nrow(NoNA.xyz),]
  maxrow <- which.max(data.warm$y)
  max <- data.warm$z[maxrow]
  data.warm <- data.warm[data.warm$z < max+1,]
  
  plot(data.cool$x, data.cool$y)
  points(data.cool$x, data.cool$y, pch=16, col="cadetblue4")
  points(data.warm$x, data.warm$y, pch=16, col="chocolate2")

coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bElev=numeric(2))
metrics_out <- data.frame(r2=numeric(2), RMSE=numeric(2), p2=numeric(2), RMSEP=numeric(2), N_sites=numeric(2), N_obvs=numeric(2))
rownames(metrics_out) <- c("Cooling", "Warming")
rownames(coeffs_out) <- c("Cooling", "Warming")


  y <- data.cool$y
  x <- data.cool$x
  z <- data.cool$z
  e <- data.cool$e
  plot(z, y)  
  plot(x, y)

  mod <- lm(y ~ x + I(x^2) + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  
  pred.y[pred.y < -0.5] = -0.5
  
  pred.new <- predict(mod, data = data.cool)
  pred.new[pred.new < -0.5] = -0.5
  data.cool$pred <- unlist(pred.new)
  
  plot(pred.y, y, main = "8-day Min Cooling Leg")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))
  gvmodel <- gvlma(mod)
  summary(gvmodel)
  plot(mod, which= 1:5)
  outlierTest(mod)
  qqPlot(mod, main="QQ Plot Cooling Leg")
  spreadLevelPlot(mod)
  plot(pred.y, mod$residuals, main="Model diagnostics Cooling Leg", xlab="Predicted", ylab="Residuals")

  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))

  library(pls)
  mod2 <- plsr(y ~ x + I(x^2) + e, validation = "LOO")
  p2 <- R2(mod2)
  detach("package:pls", unload=TRUE)


  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "Cooling"
  pred.out[,5] <- yearPath1
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_Min_", basin, "_", yrPath1, "_", yrPath2, "_sp_fall.csv"),sep = ",", col.names=F)  

  coeffs_out[1,1] <- coeffs[1,1]
  coeffs_out[1,2] <- coeffs[2,1]
  coeffs_out[1,3] <- coeffs[3,1]
  coeffs_out[1,4] <- coeffs[4,1]
  
  
  metrics_out[1,1] <- sum_mod$adj.r.squared
  metrics_out[1,2] <- post_mod$sigma
  metrics_out[1,3] <- p2$val[5]
  metrics_out[1,4] <- RMSEP
  metrics_out[1,5] <- length(unique(data.cool$SiteName))
  metrics_out[1,6] <- length(y)


  y <- data.warm$y
  x <- data.warm$x
  z <- data.warm$z
  e <- data.warm$e
  plot(z, y)  
  plot(x, y)

  mod <- lm(y ~ x + I(x^2) + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5
  
  pred.new <- predict(mod, data=data.warm)
  pred.new[pred.new < -0.5] = -0.5
  data.warm$pred <- unlist(pred.new)
  
  plot(pred.y, y, main = "8-day Min Warming Leg")
  abline(0,1)
  post_mod <- summary(lm(y ~ pred.y))

  plot(mod, which = 1:5)
  gvmodel <- gvlma(mod)
  summary(gvmodel)
  outlierTest(mod)
  qqPlot(mod, main="QQ Plot Warming Leg")
  spreadLevelPlot(mod)
  plot(pred.y, mod$residuals, main="Model diagnostics Warming Leg", xlab="Predicted", ylab="Residuals")

  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))

  library(pls)
  mod2 <- plsr(y ~ x + I(x^2) + e, validation = "LOO")
  p2 <- R2(mod2)
  detach("package:pls", unload=TRUE)


  pred.out <- matrix(nrow=length(pred.y), ncol=5)
  pred.out[,1] <- y
  pred.out[,2] <- pred.y
  pred.out[,3] <- z
  pred.out[,4] <- "Warming"
  pred.out[,5] <- yearPath2

  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_Min_", basin, "_", yrPath1, "_", yrPath2, "_sp_fall.csv"),sep = ",", col.names=F) 


  coeffs_out[2,1] <- coeffs[1,1]
  coeffs_out[2,2] <- coeffs[2,1]
  coeffs_out[2,3] <- coeffs[3,1]
  coeffs_out[2,4] <- coeffs[4,1]
  

  metrics_out[2,1] <- sum_mod$adj.r.squared
  metrics_out[2,2] <- post_mod$sigma
  metrics_out[2,3] <- p2$val[5]
  metrics_out[2,4] <- RMSEP
  metrics_out[2,5] <- length(unique(data.cool$SiteName))
  metrics_out[2,6] <- length(y)

write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", yrPath1, "_", yrPath2, "_mod_coeffs_Min.csv"), sep = ",", col.names=NA)

write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", yrPath1, "_", yrPath2, "_mod_metrics_Min.csv"), sep = ",", col.names=NA)

  pred.y <- read.csv(paste0("jk_pred_v_y_Min_", basin, "_", yrPath1, "_", yrPath2, "_sp_fall.csv"), stringsAsFactors = FALSE)
  colnames(pred.y) <- c("Y", "PredY", "JulDay", "Season", "Year")
  
  plot(pred.y$PredY, pred.y$Y, pch=16, col="blue", main="Min 8-day stream temp 2014-2015", xlab="Predicted", ylab="Observed")
  abline(0,1)
  abline(lm(pred.y$Y~ pred.y$PredY), col="blue")
  fit <- lm(pred.y$Y~ pred.y$PredY)
  plot(fit)
  summary(fit)

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day Min temp estimates
########################################################################################################

 
  setwd(paste0(mainPath, "/", subDir, "/", longBasin))
  LST.in <- read.dbf(paste0("LST_", yrPath1, "_", yrPath2, "_", basin, "_RCA.dbf"))
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  
  newnamesnum <- as.numeric(colnames(LST.in)[1:46])
  elev.in <- read.csv("Ent_rca_elev.csv", stringsAsFactors=FALSE)

  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")
  
  setwd(paste0(mainPath, subDir, longBasin, "/", modelDir, yrPath1, "_", yrPath2, "_Min/"))
  
  split.col <- which(colnames(LST.elev) == split)
  coeffs.in <- read.csv(paste0("All_data_", yrPath1, "_", yrPath2, "_mod_coeffs_Min.csv"), stringsAsFactors=FALSE)
  LogPred.out <- LST.elev[,c(1, which(colnames(LST.elev) == min): which(colnames(LST.elev) == max), 48)]

  LogPred.out[,2:(dim(LogPred.out)[2]-1)] <- 0
  rcas <- unique(elev.in$rca_id)
  
  
  for (i in 1:length(rcas))  
      {
        x <- unlist(LST.elev[i,])
        
        for (l in 2:split.col)
          {x[l] <- x[l] * coeffs.in$bLST[1] + x["Elev"] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
           }
        
        for (l in (split.col+1):dim(LogPred.out)[2]-1)     
          {x[l] <- x[l] * coeffs.in$bLST[2] + x["Elev"] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
           }
          
        fill <- rollmean(x[2:dim(LogPred.out)[2]-1], 4, align = "left")
        x[(split.col-4):(split.col+4)] <- fill[(split.col-4):(split.col+4)]
        
      LogPred.out[i,2:dim(LogPred.out)[2]-1] <- x [2:dim(LogPred.out)[2]-1] 
    }

  LogPred.out <- as.data.frame(LogPred.out)
  

plot(2:(dim(LogPred.out)[2]-1), LogPred.out[LogPred.out$RCAID == 12,2:(dim(LogPred.out)[2]-1)])
points(2:(dim(LogPred.out)[2]-1), LogPred.out[LogPred.out$RCAID == 225,2:(dim(LogPred.out)[2]-1)], pch=16, col="blue")
points(2:(dim(LogPred.out)[2]-1), LogPred.out[LogPred.out$RCAID == 371,2:(dim(LogPred.out)[2]-1)], pch=16, col="green")
points(2:(dim(LogPred.out)[2]-1), LogPred.out[LogPred.out$RCAID == 299,2:(dim(LogPred.out)[2]-1)], pch=16, col="red")
points(2:(dim(LogPred.out)[2]-1), LogPred.out[LogPred.out$RCAID == 173,2:(dim(LogPred.out)[2]-1)], pch=16, col="lightblue")
points(2:(dim(LogPred.out)[2]-1), LogPred.out[LogPred.out$RCAID == 307,2:(dim(LogPred.out)[2]-1)], pch=16, col="purple")

  LogPred.out[LogPred.out< -0.5] = -0.5
  
  newnamesnum2 <- as.numeric(colnames(LogPred.out)[2:(dim(LogPred.out)[2]-1)])
  names.out <- sprintf("Tmin_%05d", newnamesnum2)
  colnames(LogPred.out)[2:(dim(LogPred.out)[2]-1)] <- names.out
  LogPred.out$Basin_RCA <- paste0(basin, "_", LogPred.out$RCAID)

write.dbf(LogPred.out, file = paste0("predt_", yrPath1, "_", yrPath2, "_", basin, "_8D_Min.dbf"))

##########################
#This parts formats the error by day/site info
###########################

  names.out <- unique(NoNA.xyz$z)
  error.pts <- rbind(data.cool, data.warm)
  SiteID <- unique(error.pts$SiteName)
  SiteID <- as.matrix(SiteID)
  Error.pts.out <- matrix(nrow = length(SiteID), ncol = (length(names.out)+1))

  errorName <- paste0("JulDay_", names.out)
  
  Error.pts.out[1:length(SiteID), 1] <- unlist(SiteID)[1:length(SiteID)]
  Error.pts.out <- as.data.frame(Error.pts.out, stringsAsFactors = FALSE)
  colnames(Error.pts.out)[2:(length(names.out)+1)] <- names.out
  colnames(Error.pts.out)[1] <- "SiteName"
  

  for (i in SiteID) 
    { 
      error.site <- error.pts[error.pts$SiteName  == i,]
      error.site$error <- error.site[,'y']-error.site[,'pred']
      
      
      error <- matrix(ncol=1, nrow=length(names.out))
      error[,1] <- names.out
      error <- data.frame(error)
      colnames(error) <- c("JulDay")
      
      error.site.fill <- merge(error, error.site, by.x = "JulDay", by.y = "z", all.x=TRUE, all.y = FALSE)
      Error.pts.out[Error.pts.out$SiteName==i,2:(length(names.out)+1)] <- as.numeric(unlist(error.site.fill$error))
    }

  Error.pts.out[,2:(length(names.out)+1)] <- sapply(Error.pts.out[,2:(length(names.out)+1)], as.numeric)

  Error.pts.out[,2:(length(names.out)+1)] <- round(Error.pts.out[,2:(length(names.out)+1)], digits=3)

write.dbf(Error.pts.out, file = paste0("Error", "_", basin, "_", yrPath1, "_", yrPath2, "_8D_Min.dbf")) 
write.csv(Error.pts.out, file = paste0("Error", "_", basin, "_", yrPath1, "_", yrPath2, "_8D_Min.csv"))

###############################################
# shapefile output
###############################################
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
yrPath1 <- "14"
yrPath2 <- "15"
dataPath <- "D:/OneDrive/work/research/CHaMP/GIS/coverages/"
mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
seis <- rev(seis)

  projname <- paste0(longBasin, "_net_rca")

  setwd(paste0(dataPath, longBasin))

  proj_layer <- readOGR(dsn=".", layer = projname)

  setwd(paste0(mainPath, "/Min_temp_models/", yrPath1, "_", yrPath2, "/"))
  
  netname <- paste0(basin, "_", yrPath1, "_", yrPath2, "_8D_winter_min")
  network <- readOGR(dsn=".", layer = netname)
  network <- spTransform(network, proj4string(proj_layer))

# #### plot it to make sure it looks right ######### 


  fix3 <- classIntervals(network@data[,10], n = 11, style = "fixed",fixedBreaks=c(-1,0,1,2,4,6,8,10,12,14,16))
  fix3.colors <- findColours(fix3,pal=seis)
  plot(network, col=fix3.colors, bg="black", fg="white")

################################################
# animation
################################################

library(rgdal)
library(RColorBrewer)
library(classInt)

setwd(paste0(mainPath, "/Min_temp_models/", yrPath1, "_", yrPath2, "/"))
  
modelPath <- paste0(mainPath, "Min_temp_models/")

shpPath <- paste0(modelPath, "graphics/")

  netname <- paste0(basin, "_", yrPath1, "_", yrPath2, "_8D_winter_min")
  ptsname <- paste0(basin, "_", yrPath1, "_", yrPath2, "_8D_winter_min_error")

  error_pts <- readOGR(dsn=".", ptsname)
  error_pts <- spTransform(error_pts, proj4string(proj_layer))
  
  network <- readOGR(dsn=".", layer = netname)


setwd(shpPath)

  seis = c("#AA0000", "#D00000", "#F70000", "#FF1D00", "#FF4400", "#FF6A00", "#FF9000", "#FFB700", "#FFDD00", "#FFE200", "#BDFF0C", "#73FF1A", "#3FFA36", "#16F45A", "#00D08B", "#0087CD", "#0048FA", "#0024E3")
  seis <- rev(seis)

  names.out <- colnames(network@data[3:39])
  namesnum <- as.numeric(gsub("Tmin_", "", colnames(network@data[3:39])))
  means <- colMeans(network@data[3:39])
  SDs <- colStdevs(network@data[3:39])
  yplus <- means + SDs
  yminus <- means - SDs
  df <- data.frame(means=means, SDs=SDs, names=namesnum)
  sequ <- c(1:37)
  namer <- sprintf('%03d', sequ)
  fix4 <- classIntervals(means, n = 10, style = "fixed",fixedBreaks=c(-1,0,1,2,4,6,8,10,12,14,16))
  fix4.colors <- findColours(fix4,pal=seis)

  for (i in 3:39)
    {
      
      namey <- gsub("Tmin_", "", colnames(network@data)[i])
      
      filename <- paste0(shpPath, namer[i-2], ".png", sep="")
      png(filename=filename, res = 300, width = 1500, height = 1500, units = "px", bg="black")
      
      fix3 <- classIntervals(network@data[,i], n = 11, style = "fixed",fixedBreaks=c(-1,0,1,2,4,6,8,10,12,14,20))
      fix3.colors <- findColours(fix3,pal=seis)
      
      cexEr <- ifelse(abs(error_pts@data[,i-1]) <= 1, 0.5,
                      ifelse(abs(error_pts@data[,i-1])>1, 0.75,
                             ifelse(abs(error_pts@data[,i-1])>2, 1.0,
                                    ifelse(abs(error_pts@data[,i-1])>3, 1.25, NA))))
      
      plot(network, col=fix3.colors, bg="black", fg="white")
      points(error_pts, pch=16, col="gray40", cex=cexEr)
      
      legend("right", fill = attr(fix3.colors, "palette"), title="8-day Min (°C)", legend = c("-1-0","0-1","1-2","2-3","3-4","4-6","6-8","8-10","10-12","12+"), bty = "n", cex=.5, inset=c(.1,0), text.col="white");
      legend(x=grconvertX(c(0.05, 0.25), from='npc'), 
             y=grconvertY(c(0.3, 0.5), from='npc'),, title="Model error (°C)", legend = c("0-1","1-2","2-3","3+"), bty = "n", pch=16, pt.cex=c(0.5, 0.75, 1.0, 1.5), col = "gray40", cex=.5, text.col="white");
      
      
      title("Entiat 8-day min 2014-2015 (°C)", sub = paste0("Julian Day ", namey), line=-0.9, adj=.80, col.main="white", col.sub="white", outer=FALSE, cex.main=0.5, cex.sub=0.5)
      tmp2 <- subplot(
        plot(sequ[1:(i-2)], means[1:(i-2)], col=fix4.colors, pch=16, bty="n", xlim=c(1,44), ylim=c(-1,18), cex.main=.8, main="Basin mean", adj=0, xaxt='n', xlab='',ylab='', col.lab="white", cex.axis=0.5, cex.lab = 0.25, col.axis="white", col.main = "white", bg="black"), 
        x=grconvertX(c(0.1,0.45), from='npc'), 
        y=grconvertY(c(0.05, 0.20), from='npc'),
        size=c(1,1.5), vadj=0.5, hadj=0.5, 
        pars=list( mar=c(0,0,0,0)+0.1, cex=0.5))
      
      
      dev.off()
    }
  

setwd(shpPath) 
system('"C:/Program Files/ImageMagick-7.0.1-Q16/convert.exe" -delay 20 -morph 3 *.png Ent_14_15_min.mpeg')
