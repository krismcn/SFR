############################################################################################################
# This set of R scripts combines a series of scripts that make up the work flow to predict and validate 
# stream temperature for a year using MODIS 1km LST data
# The initial input is an LST_YY.dbf which is a table (columns = days, rows = sites) with the grid value for each grid cell in the spatial extent
# Should do a search-and-replace for the output folder (usually dated), and the year being processed, in both a YYYY and a _YY format.

# Edited Aug 2014 to add the PRESS stastic output
# Edited Aug 2015 to refine code flow and documentation

          
##############################################################################################
# This section reads in the 1km LST data for a year and fills any gaps across the year at each pixel 
# using a linear interpolation.
##############################################################################################
setwd("D:/Dropbox/work/research/CHaMP/CHaMP_data/Asotin")

library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

	LST.in <- read.dbf("Aso_LST_2013.dbf")
	GrPolID <- LST.in[,1]
	LST.in <- LST.in[,3:48]
	LST.in[LST.in<0.1] = NA
	
	LST.names<-colnames(LST.in)
	LST.names <- order(LST.names)
	LST.in<- LST.in[,(LST.names)]
	
	tLST.in <- t(LST.in)
	tLST.out <- interpNA(tLST.in, method = "linear")
	CLST <- tLST.out*0.02-273.15
	tCLST <- t(CLST)
	CLST.out <- as.data.frame(tCLST)
	plot(1:46, CLST.out[100,1:46])

	myFunc <- function(x)	#in case the first day is missing data
	            {
	            if (is.na(x[1]))
	            {
	            x[1] <- x[2]
	            }
	            x
	            }
            
	LST.filled <- apply(CLST.out, 1, function(x) myFunc(x))

	myFunc2 <- function(x)	#in case the second to last day is missing data
	{
	  if (is.na(x[45]))
	  {
	    x[45] <- x[44]
	  }
	  x
	}

	tLST.filled <- t(LST.filled)
	LST.filled <- apply(tLST.filled, 1, function(x) myFunc2(x))

	myFunc3 <- function(x)	#in case the last day is missing data
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
	LST.out[,1:46] <- round(LST.out[,1:46], digits = 2)
	
	write.dbf(LST.out, file = "LST13_Aso_interp.dbf")

####################################################################
#This part calculates the zonal mean of the LST for each RCA
####################################################################

	out_Preds <- read.dbf("LST13_Aso_interp.dbf") 

	weights <- read.dbf("Aso_rca_area_wgts.dbf")	#these are the area weights of each polygon from the intersected 1km poly and rca shapefiles
	weights[weights$area_wgts<0.01, "area_wgts"] = 0.01
	rcas <- unique(unlist(weights$RCAID))
	rca_zonal <- matrix(ncol=46, nrow=length(rcas))

	 l <- 1
	  for(i in rcas)	
	      {
	      pixels <- weights[weights$RCAID == i, "GRIDCODE"]
	      wgts <- weights[weights$RCAID == i, "area_wgts"]
	      
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
	
	
	write.dbf(rca_zonal, file = "LST14_Aso_RCA.dbf")



#################################################################
#Logger prediction modeling part
#Generates models for all model structures and selects "best" model using the PRESS statistic
#################################################################

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/Asotin/"
yearPath <- "2013"
setwd(paste0(mainPath))

  ID.in <- read.csv("Aso_sites_GRIDCODE.csv")
  
  Log.in <- read.csv("Aso_8Day_Mn_2013.csv")
  
   
  LST.in <- read.dbf("LST13_Aso_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- substring (colnames(LST.in)[1:46],3,5)
  newnamesnum <- as.numeric(newnamesnum)
  colnames(LST.in)[1:46] <- newnamesnum

  SiteID <- unique(Log.in$OrigID)
  SiteID <- as.matrix(SiteID)

#_____________________________________________________________________
# Use this secion if daily logger data needs to be summarized to 8day first 
#_____________________________________________________________________
  dates <- matrix(nrow=365, ncol=1)
  dates[,1] <- 1:365
  colnames(dates) <- c("Days")

 
  Log.8Day.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.in[Log.in$OrigID == i,]
      full.year <- merge(dates, Log.site, by.x = "Days", by.y = "JulDay", all.x = TRUE)
      eightday <- rollapply(full.year$DailyMaximum, width = 8, FUN = mean, na.rm = T, fill = NA)
      eightday <- as.matrix(eightday)
      full.year$Mx8D <- eightday
      full.year$Mx8D[361] <- mean(full.year$DailyMaximum[361:365])
      Log.8Day.out <- rbind(Log.8Day.out, full.year)
      
    }

  Log.in <- Log.8Day.out
#_________________________________________________________________________

setwd(paste0(mainPath, yearPath))

  
  LST.Log.out <- data.frame (mup = NULL)
  
    for (i in SiteID) 
      { 
        Log.site <- Log.in[Log.in$SiteName == i,]
        Log.site <- as.data.frame(Log.site)
        
        
        RCAID <- ID.in$rca_id[ID.in$SiteName == i]
        Elev <- ID.in$ELEV[ID.in$SiteName == i]
        
        LST.site <- matrix(ncol=3, nrow=46)
        LST.site[,1] <- as.numeric(unlist(colnames(LST.in)[1:46]))
        LST.site[,2] <- unlist(LST.in[RCAID,1:46])
        LST.site <- data.frame(LST.site)
        colnames(LST.site) <- c("JulDay", "LST", "Elev")
        LST.site[3] <- Elev
        LST.Log.site <- merge(LST.site, Log.site, by.x = "JulDay", by.y = "Days", all.x=TRUE, all.y = FALSE)
        
        LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
      }
  
  
  
  ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
  NoNA.xyz <- LST.Log.out[ind,]
  NoNA.xyz <- NoNA.xyz[,c(15, 2, 1, 3, 7)]
  colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
  plot(NoNA.xyz$x, NoNA.xyz$y)
  
  write.csv(x=NoNA.xyz, file="Aso_2013_8Day_model_data.csv", row.names = FALSE)
  
  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  
  maxrow <- which.max(NoNA.xyz$y)
  data.sp <- NoNA.xyz[1:maxrow,]
  data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]


#################################
# summer data only models
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
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Aso_2013_summer_only.csv",sep = ",", col.names=F)  

  plot(pred.out[,1], pred.out[,2])
  abline(0,1)

################################
# full year data model
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
  pred.out[,5] <- "2013"
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Aso_2013_full_year.csv",sep = ",", col.names=F)  

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
  
  write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_Aso_2013_sp_fall.csv",sep = ",", col.names=F)  
  
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
# This part applies the model coefficients to the LST to generate predicted stream temps
# 
########################################################################################################

#first section has code for basin-wide model application

mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/Asotin/"
yearPath <- "2013"
setwd(paste0(mainPath))

  
  LST.in <- read.dbf("LST13_Aso_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- substring (colnames(LST.in)[1:46],3,5)
  newnamesnum <- as.numeric(newnamesnum)
  colnames(LST.in)[1:46] <- newnamesnum
  

    
  elev.in <- read.dbf("Aso_rca_elev.dbf")
  elev.in$RCAID <- as.numeric(gsub("Asotin_", "", elev.in$BASINRCA))
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCAID")
  LST.elev <- LST.elev[,c(1:47,52)]
  colnames(LST.elev)[48] <- "elev"
 
setwd(paste0(mainPath, yearPath))

  coeffs.in <- read.csv("All_data_2013_mod_coeffs_Mx.csv")
  LogPred.out <- LST.elev
  LogPred.out[,2:48] <- 0
  rcas <- unique(elev.in$RCAID)


  for (i in rcas)
    {
      x <- unlist(LST.elev[LST.elev$RCAID == i,2:48])
      maxrow <- as.numeric(which.max(x[1:12])) #either specify or let be dynamic
      day <- as.numeric(colnames(LST.elev)[maxrow])
      
      j <- 1
      
      for (l in 1:maxrow)
        {x[l] <- x[l] * coeffs.in$sp_bLST[5] + j * coeffs.in$sp_bJul[5] +  x[47] * coeffs.in$sp_bElev[5] + x[l]^2 * coeffs.in$sp_bLST2[5] + coeffs.in$sp_Int[5]
         j <- j + 8}
      k <- day
      for (l in midrow:46)     
        {x[l] <- x[l] * coeffs.in$fall_bLST[5] + k * coeffs.in$fall_bJul[5] + x[47] * coeffs.in$fall_bElev[5] + x[l]^2 * coeffs.in$fall_bLST2[5] + coeffs.in$fall_Int[5]
         k <- k + 8}
      LogPred.out[i,2:48] <- x  
    }


  LogPred.out <- as.data.frame(LogPred.out)
  LogPred.out[LogPred.out< -0.5] = -0.5
  
  names.out <- sprintf("Tmn_13_%03d", newnamesnum)
  colnames(LogPred.out)[2:47] <- names.out
  LogPred.out$Basin <- "Asotin"
  colnames(LogPred.out)[1] <- "RCAID"
  LogPred.out$BasinRCA <- paste0(LogPred.out$Basin, "_", LogPred.out$RCAID) 
  LogPred.out$vModel <- "11Aug2015"
write.dbf(LogPred.out, file = "predt2013_Aso_8D_Mn.dbf")

# DONE




######################################################################################
# this section has a function for more complicated basin partitioning model application
# in this case - partitioning by HUC
######################################################################################

HUCs.in <- read.csv("D:/Dropbox/work/research/Steelhead/spatial_data/coverages/RCA_HUCs.csv")
LST.HUCs <- merge(LST.in, HUCs.in, by.x = "OBJECTID", by.y = "OBJECTID", all.x = F, all.y = T)
LST.HUCs.out <- LST.HUCs
LST.HUCs.out[,] <- 0


myFunc <- function(x) #HUC based model application
{
  maxrow <- as.numeric(which.max(x[2:47]))
  midrow <- maxrow + 1
  
  if (x[48] == 1)
  {
    for (i in 2:maxrow)
    {x[i] <- x[i] * coeffs.in$sp_noJul_bLST[1] + x[i]^2 * coeffs.in$sp_noJul_bLST2[1] + coeffs.in$sp_noJul_Int[1]}
    for (i in midrow:47)     
    {x[i] <- x[i] * coeffs.in$fall_noJul_bLST[1] + x[i]^2 * coeffs.in$fall_noJul_bLST2[1]  + coeffs.in$fall_noJul_Int[1]} 
  } else
    if (x[48] == 2)
    {
      for (i in 2:maxrow)
      {x[i] <- x[i] * coeffs.in$sp_noJul_bLST[2] + x[i]^2 * coeffs.in$sp_noJul_bLST2[2]  + coeffs.in$sp_noJul_Int[2]}
      for (i in midrow:47)     
      {x[i] <- x[i] * coeffs.in$fall_noJul_bLST[2] + x[i]^2 * coeffs.in$fall_noJul_bLST2[2] + coeffs.in$fall_noJul_Int[2]} 
    } else
      if (x[48] == 3)
      {
        for (i in 2:maxrow)
        {x[i] <- x[i] * coeffs.in$sp_noJul_bLST[3] + x[i]^2 * coeffs.in$sp_noJul_bLST2[3] + coeffs.in$sp_noJul_Int[3]}
        for (i in midrow:47)     
        {x[i] <- x[i] * coeffs.in$fall_noJul_bLST[3] + x[i]^2 * coeffs.in$fall_noJul_bLST2[3] + coeffs.in$fall_noJul_Int[3]} 
      } else
        if (x[48] == 4)
        {
          for (i in 2:maxrow)
          {x[i] <- x[i] * coeffs.in$sp_noJul_bLST[4] + x[i]^2 * coeffs.in$sp_noJul_bLST2[4] + coeffs.in$sp_noJul_Int[4]}
          for (i in midrow:47)     
          {x[i] <- x[i] * coeffs.in$fall_noJul_bLST[4] + x[i]^2 * coeffs.in$fall_noJul_bLST2[4] + coeffs.in$fall_noJul_Int[4]} 
        }
  x
}

LogPred.out <- apply(LST.HUCs, 1, myFunc)
tLogPred.out <- t(LogPred.out)
tLogPred.out[tLogPred.out<0] = 0.0

tLogPred.out[,2:47] <- round(tLogPred.out[,2:47], digits = 2)
#crosswalk.in <- read.csv("D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/ObjID_crosswalk_JD_streams_alb.csv")
merged.LogPred <- merge(tLogPred.out, crosswalk.in, by.x = "OBJECTID", by.y = "OBJECTID", all.x = TRUE, all.y = TRUE)

write.dbf(merged.LogPred, file = "predt2013_JD_8D_Mn.dbf")








######################################################################################
# Older modeling script for considering different linear model structures and then 
# choosing the "best" based on the PRESS statistic before writing out coefficients
# Can be used in new watersheds if there's a need to investigate
######################################################################################

  coeffs_out <- data.frame(sp_Int=numeric(5), sp_bLST=numeric(5), sp_bLST2=numeric(5), sp_bJul=numeric(1), sp_bElev=numeric(5), fall_Int=numeric(5), fall_bLST=numeric(5), fall_bLST2=numeric(5), fall_bJul=numeric(5), fall_bElev=numeric(5))
  metrics_out <- data.frame(PRESS_sp=numeric(5), p2_sp=numeric(5), r2_sp=numeric(5), RMSEP_sp=numeric(5), RMSE_sp=numeric(5), AIC_sp=numeric(5), PRESS_fall=numeric(5), p2_fall=numeric(5), r2_fall=numeric(5), RMSEP_fall=numeric(5), RMSE_fall=numeric(5), AIC_fall=numeric(5))
  rownames(metrics_out) <- c("noJul", "Jul", "Poly", "Elev", "All")
  rownames(coeffs_out) <- c("noJul", "Jul", "Poly", "Elev", "All")

  y.spring <- NoNA.xyz$Avg8Day[1:maxrow]
  x.spring <- NoNA.xyz$LST[1:maxrow]
  z.spring <- NoNA.xyz$JulDay[1:maxrow]
  e.spring <- NoNA.xyz$Elev[1:maxrow]


  noJul_mod_sp <- lm(y.spring ~ x.spring)
  noJul_summ_sp <- summary(lm(y.spring ~ x.spring))
  noJul_plsr_mod_sp <- plsr(y.spring ~ x.spring, validation = "LOO")
  noJul_RMSEP_sp <- RMSEP(noJul_plsr_mod_sp)
  sp.noJul.coeffs <- coefficients(noJul_summ_sp)
  
  Jul_mod_sp <- lm(y.spring ~ x.spring + z.spring)
  Jul_summ_sp <- summary(lm(y.spring ~ x.spring + z.spring))
  Jul_plsr_mod_sp <- plsr(y.spring ~ x.spring + z.spring, validation = "LOO")
  Jul_RMSEP_sp <- RMSEP(Jul_plsr_mod_sp)
  sp.Jul.coeffs <- coefficients(Jul_summ_sp)
  
  poly_mod_sp <- lm(y.spring ~ x.spring + I(x.spring^2))
  poly_summ_sp <- summary(lm(y.spring ~ x.spring + I(x.spring^2)))
  poly_plsr_mod_sp <- plsr(y.spring ~ x.spring + I(x.spring^2) + z.spring , validation = "LOO")
  poly_RMSEP_sp <- RMSEP(poly_plsr_mod_sp)
  sp.poly.coeffs <- coefficients(poly_summ_sp)
  
  elev_mod_sp <- lm(y.spring ~ x.spring + z.spring + e.spring)
  elev_summ_sp <- summary(lm(y.spring ~ x.spring + z.spring + e.spring))     
  elev_plsr_mod_sp <- plsr(y.spring ~ x.spring + z.spring + e.spring, validation = "LOO")
  elev_RMSEP_sp <- RMSEP(elev_plsr_mod_sp)
  sp.elev.coeffs <- coefficients(elev_summ_sp)
  
  all_mod_sp <- lm(y.spring ~ x.spring + I(x.spring^2) + z.spring + e.spring)
  all_summ_sp <- summary(lm(y.spring ~ x.spring + I(x.spring^2) + z.spring + e.spring))   
  all_plsr_mod_sp <- plsr(y.spring ~ x.spring + I(x.spring^2) + z.spring + e.spring, validation = "LOO")
  all_RMSEP_sp <- RMSEP(all_plsr_mod_sp)    
  sp.all.coeffs <- coefficients(all_summ_sp)
  
  pressstat_Jul <- PRESS(Jul_mod_sp, verbose = "FALSE")
  pressstat_noJul <- PRESS(noJul_mod_sp, verbose = "FALSE")
  pressstat_poly <- PRESS(poly_mod_sp, verbose = "FALSE")
  pressstat_elev <- PRESS(elev_mod_sp, verbose = "FALSE")
  pressstat_all <- PRESS(all_mod_sp, verbose = "FALSE")
  
  pressstat_Jul$stat
  pressstat_noJul$stat
  pressstat_poly$stat
  pressstat_elev$stat
  pressstat_all$stat
  
  AIC(noJul_mod_sp)
  AIC(Jul_mod_sp)
  AIC(poly_mod_sp)
  AIC(elev_mod_sp)
  AIC(all_mod_sp)
  
  
  
  metrics_out[1,1] <- pressstat_noJul$stat
  metrics_out[1,2] <- pressstat_noJul$P.square
  metrics_out[1,3] <- noJul_summ_sp$adj.r.squared     
  metrics_out[1,4] <- noJul_RMSEP_sp[[1]][[4]]
  metrics_out[1,5] <- noJul_summ_sp$sigma
  metrics_out[1,6] <- AIC(noJul_mod_sp)
  metrics_out[2,1] <- pressstat_Jul$stat
  metrics_out[2,2] <- pressstat_Jul$P.square
  metrics_out[2,3] <- Jul_summ_sp$adj.r.squared     
  metrics_out[2,4] <- Jul_RMSEP_sp[[1]][[4]]
  metrics_out[2,5] <- Jul_summ_sp$sigma
  metrics_out[2,6] <- AIC(Jul_mod_sp)
  metrics_out[3,1] <- pressstat_poly$stat
  metrics_out[3,2] <- pressstat_poly$P.square
  metrics_out[3,3] <- poly_summ_sp$adj.r.squared     
  metrics_out[3,4] <- poly_RMSEP_sp[[1]][[4]]
  metrics_out[3,5] <- poly_summ_sp$sigma
  metrics_out[3,6] <- AIC(poly_mod_sp)
  metrics_out[4,1] <- pressstat_elev$stat
  metrics_out[4,2] <- pressstat_elev$P.square
  metrics_out[4,3] <- elev_summ_sp$adj.r.squared     
  metrics_out[4,4] <- elev_RMSEP_sp[[1]][[4]]
  metrics_out[4,5] <- elev_summ_sp$sigma
  metrics_out[4,6] <- AIC(elev_mod_sp)
  metrics_out[5,1] <- pressstat_all$stat
  metrics_out[5,2] <- pressstat_all$P.square
  metrics_out[5,3] <- all_summ_sp$adj.r.squared     
  metrics_out[5,4] <- all_RMSEP_sp[[1]][[4]]
  metrics_out[5,5] <- all_summ_sp$sigma
  metrics_out[5,6] <- AIC(all_mod_sp)
  
  
  if (pressstat_Jul$stat < pressstat_noJul$stat & pressstat_Jul$stat < pressstat_poly$stat & pressstat_Jul$stat < pressstat_elev$stat & pressstat_Jul$stat < pressstat_all$stat)
  {pred.y.spring<- predict(noJul_mod_sp); sp_mod <- noJul_mod_sp}
  pressstat_Jul$stat < pressstat_noJul$stat & pressstat_Jul$stat < pressstat_poly$stat & pressstat_Jul$stat < pressstat_elev$stat & pressstat_Jul$stat < pressstat_all$stat
  
  if (pressstat_noJul$stat < pressstat_Jul$stat & pressstat_noJul$stat < pressstat_poly$stat & pressstat_noJul$stat < pressstat_elev$stat & pressstat_noJul$stat < pressstat_all$stat)
  {pred.y.spring<- predict(Jul_mod_sp); sp_mod <- Jul_mod_sp}
  pressstat_noJul$stat < pressstat_Jul$stat & pressstat_noJul$stat < pressstat_poly$stat & pressstat_noJul$stat < pressstat_elev$stat & pressstat_noJul$stat < pressstat_all$stat      
  
  if (pressstat_poly$stat < pressstat_noJul$stat & pressstat_poly$stat < pressstat_Jul$stat & pressstat_poly$stat < pressstat_elev$stat & pressstat_poly$stat < pressstat_all$stat)
  {pred.y.spring<- predict(poly_mod_sp); sp_mod <- poly_mod_sp}
  pressstat_poly$stat < pressstat_noJul$stat & pressstat_poly$stat < pressstat_Jul$stat & pressstat_poly$stat < pressstat_elev$stat & pressstat_poly$stat < pressstat_all$stat
  
  if (pressstat_elev$stat < pressstat_noJul$stat & pressstat_elev$stat < pressstat_Jul$stat & pressstat_elev$stat < pressstat_poly$stat & pressstat_elev$stat < pressstat_all$stat)
  {pred.y.spring<- predict(elev_mod_sp); sp_mod <- elev_mod_sp}
  pressstat_elev$stat < pressstat_noJul$stat & pressstat_elev$stat < pressstat_Jul$stat & pressstat_elev$stat < pressstat_poly$stat & pressstat_elev$stat < pressstat_all$stat      
  
  if (pressstat_all$stat < pressstat_noJul$stat & pressstat_all$stat < pressstat_Jul$stat & pressstat_all$stat < pressstat_poly$stat & pressstat_all$stat < pressstat_elev$stat)
  {pred.y.spring<- predict(all_mod_sp); sp_mod <- all_mod_sp}
  pressstat_all$stat < pressstat_noJul$stat & pressstat_all$stat < pressstat_Jul$stat & pressstat_all$stat < pressstat_poly$stat & pressstat_all$stat < pressstat_elev$stat      
  
  pred.y.spring[pred.y.spring<0.5] = 0
  pred.out <- matrix(nrow=length(y.spring), ncol=4)
  colnames(pred.out) <- c("Y_sp", "PredY_sp", "SiteName", "Year")
  pred.out[,1] <- y.spring
  pred.out[,2] <- pred.y.spring
  pred.out[,3] <- as.character(NoNA.xyz$SiteName[1:maxrow])
  pred.out[,4] <- 2013
  plot(pred.out[,1], pred.out[,2])
  
  setwd(pathname)
  
  write.table(x=pred.out, append=F,row.names=F, file = "2013_PredY_Sp.csv", sep = ",", col.names=T)
  
  
  
  y.fall <- NoNA.xyz$Avg8Day[maxrow:nrow(NoNA.xyz)]
  x.fall <- NoNA.xyz$LST[maxrow:nrow(NoNA.xyz)]
  z.fall <- NoNA.xyz$JulDay[maxrow:nrow(NoNA.xyz)]
  e.fall <- NoNA.xyz$Elev[maxrow:nrow(NoNA.xyz)]
  plot(x.fall, y.fall)
  
  noJul_mod_fall <- lm(y.fall ~ x.fall)
  noJul_summ_fall <- summary(lm(y.fall ~ x.fall))
  noJul_plsr_mod_fall <- plsr(y.fall ~ x.fall, validation = "LOO")
  noJul_RMSEP_fall <- RMSEP(noJul_plsr_mod_fall)
  fall.noJul.coeffs <- coefficients(noJul_summ_fall)
  
  Jul_mod_fall <- lm(y.fall ~ x.fall + z.fall)
  Jul_summ_fall <- summary(lm(y.fall ~ x.fall + z.fall))
  Jul_plsr_mod_fall <- plsr(y.fall ~ x.fall + z.fall, validation = "LOO")
  Jul_RMSEP_fall <- RMSEP(Jul_plsr_mod_fall)
  fall.Jul.coeffs <- coefficients(Jul_summ_fall)
  
  poly_mod_fall <- lm(y.fall ~ x.fall + I(x.fall^2))
  poly_summ_fall <- summary(lm(y.fall ~ x.fall + I(x.fall^2)))
  poly_plsr_mod_fall <- plsr(y.fall ~ x.fall + I(x.fall^2) + z.fall , validation = "LOO")
  poly_RMSEP_fall <- RMSEP(poly_plsr_mod_fall)
  fall.poly.coeffs <- coefficients(poly_summ_fall)
  
  elev_mod_fall <- lm(y.fall ~ x.fall + z.fall + e.fall)
  elev_summ_fall <- summary(lm(y.fall ~ x.fall + z.fall + e.fall))     
  elev_plsr_mod_fall <- plsr(y.fall ~ x.fall + z.fall + e.fall, validation = "LOO")
  elev_RMSEP_fall <- RMSEP(elev_plsr_mod_fall)
  fall.elev.coeffs <- coefficients(elev_summ_fall)
  
  all_mod_fall <- lm(y.fall ~ x.fall + I(x.fall^2) + z.fall + e.fall)
  all_summ_fall <- summary(lm(y.fall ~ x.fall + I(x.fall^2) + z.fall + e.fall))   
  all_plsr_mod_fall <- plsr(y.fall ~ x.fall + I(x.fall^2) + z.fall + e.fall, validation = "LOO")
  all_RMSEP_fall <- RMSEP(all_plsr_mod_fall)    
  fall.all.coeffs <- coefficients(all_summ_fall)
  
  pressstat_Jul <- PRESS(Jul_mod_fall, verbose = "FALSE")
  pressstat_noJul <- PRESS(noJul_mod_fall, verbose = "FALSE")
  pressstat_poly <- PRESS(poly_mod_fall, verbose = "FALSE")
  pressstat_elev <- PRESS(elev_mod_fall, verbose = "FALSE")
  pressstat_all <- PRESS(all_mod_fall, verbose = "FALSE")
  
  pressstat_Jul$stat
  pressstat_noJul$stat
  pressstat_poly$stat
  pressstat_elev$stat
  pressstat_all$stat
  
  AIC(noJul_mod_fall)
  AIC(Jul_mod_fall)
  AIC(poly_mod_fall)
  AIC(elev_mod_fall)
  AIC(all_mod_fall)
  
  
  metrics_out[1,7] <- pressstat_noJul$stat
  metrics_out[1,8] <- pressstat_noJul$P.square
  metrics_out[1,9] <- noJul_summ_fall$adj.r.squared     
  metrics_out[1,10] <- noJul_RMSEP_fall[[1]][[4]]
  metrics_out[1,11] <- noJul_summ_fall$sigma
  metrics_out[1,12] <- AIC(noJul_mod_fall)
  metrics_out[2,7] <- pressstat_Jul$stat
  metrics_out[2,8] <- pressstat_Jul$P.square
  metrics_out[2,9] <- Jul_summ_fall$adj.r.squared     
  metrics_out[2,10] <- Jul_RMSEP_fall[[1]][[4]]
  metrics_out[2,11] <- Jul_summ_fall$sigma
  metrics_out[2,12] <- AIC(Jul_mod_fall)
  metrics_out[3,7] <- pressstat_poly$stat
  metrics_out[3,8] <- pressstat_poly$P.square
  metrics_out[3,9] <- poly_summ_fall$adj.r.squared     
  metrics_out[3,10] <- poly_RMSEP_fall[[1]][[4]]
  metrics_out[3,11] <- poly_summ_fall$sigma
  metrics_out[3,12] <- AIC(poly_mod_fall)
  metrics_out[4,7] <- pressstat_elev$stat
  metrics_out[4,8] <- pressstat_elev$P.square
  metrics_out[4,9] <- elev_summ_fall$adj.r.squared     
  metrics_out[4,10] <- elev_RMSEP_fall[[1]][[4]]
  metrics_out[4,11] <- elev_summ_fall$sigma
  metrics_out[4,12] <- AIC(elev_mod_fall)
  metrics_out[5,7] <- pressstat_all$stat
  metrics_out[5,8] <- pressstat_all$P.square
  metrics_out[5,9] <- all_summ_fall$adj.r.squared     
  metrics_out[5,10] <- all_RMSEP_fall[[1]][[4]]
  metrics_out[5,11] <- all_summ_fall$sigma
  metrics_out[5,12] <- AIC(all_mod_fall)
  
  if (pressstat_Jul$stat < pressstat_noJul$stat & pressstat_Jul$stat < pressstat_poly$stat & pressstat_Jul$stat < pressstat_elev$stat & pressstat_Jul$stat < pressstat_all$stat)
  {pred.y.fall<- predict(lm(y.fall ~ x.fall + z.fall)); fall_mod <- Jul_mod_fall}
  pressstat_Jul$stat < pressstat_noJul$stat & pressstat_Jul$stat < pressstat_poly$stat & pressstat_Jul$stat < pressstat_elev$stat & pressstat_Jul$stat < pressstat_all$stat
  
  if (pressstat_noJul$stat < pressstat_Jul$stat & pressstat_noJul$stat < pressstat_poly$stat & pressstat_noJul$stat < pressstat_elev$stat & pressstat_noJul$stat < pressstat_all$stat)
  {pred.y.fall<- predict(lm(y.fall ~ x.fall)); fall_mod <- noJul_mod_fall}
  pressstat_noJul$stat < pressstat_Jul$stat & pressstat_noJul$stat < pressstat_poly$stat & pressstat_noJul$stat < pressstat_elev$stat & pressstat_noJul$stat < pressstat_all$stat      
  
  if (pressstat_poly$stat < pressstat_noJul$stat & pressstat_poly$stat < pressstat_Jul$stat & pressstat_poly$stat < pressstat_elev$stat & pressstat_poly$stat < pressstat_all$stat)
  {pred.y.fall<- predict(lm(y.fall ~ x.fall + I(x.fall^2))); fall_mod <- poly_mod_fall}
  pressstat_poly$stat < pressstat_noJul$stat & pressstat_poly$stat < pressstat_Jul$stat & pressstat_poly$stat < pressstat_elev$stat & pressstat_poly$stat < pressstat_all$stat
  
  if (pressstat_elev$stat < pressstat_noJul$stat & pressstat_elev$stat < pressstat_Jul$stat & pressstat_elev$stat < pressstat_poly$stat & pressstat_elev$stat < pressstat_all$stat)
  {pred.y.fall<- predict(lm(y.fall ~ x.fall + z.fall + e.fall)); fall_mod <- elev_mod_fall}
  pressstat_elev$stat < pressstat_noJul$stat & pressstat_elev$stat < pressstat_Jul$stat & pressstat_elev$stat < pressstat_poly$stat & pressstat_elev$stat < pressstat_all$stat      
  
  if (pressstat_all$stat < pressstat_noJul$stat & pressstat_all$stat < pressstat_Jul$stat & pressstat_all$stat < pressstat_poly$stat & pressstat_all$stat < pressstat_elev$stat)
  {pred.y.fall<- predict(lm(y.fall ~ x.fall + I(x.fall^2) + z.fall + e.fall)); fall_mod <- all_mod_fall}
  pressstat_all$stat < pressstat_noJul$stat & pressstat_all$stat < pressstat_Jul$stat & pressstat_all$stat < pressstat_poly$stat & pressstat_all$stat < pressstat_elev$stat      
  
  pred.y.fall[pred.y.fall<0.5] = 0
  pred.out <- matrix(nrow=length(y.fall), ncol=4)
  colnames(pred.out) <- c("Y_fall", "PredY_fall", "SiteName", "Year")
  pred.out[,1] <- y.fall
  pred.out[,2] <- pred.y.fall
  pred.out[,3] <- as.character(NoNA.xyz$SiteName[maxrow:nrow(NoNA.xyz)])
  pred.out[,4] <- 2013
  plot(pred.out[,1], pred.out[,2])
  write.table(x=pred.out, append=F,row.names=F, file = "2013_PredY_Fall.csv", sep = ",", col.names=T)
  
  
  
  coeffs_out[1,1] <- sp.noJul.coeffs[1,1]
  coeffs_out[1,2] <- sp.noJul.coeffs[2,1]  
  coeffs_out[1,6] <- fall.noJul.coeffs[1,1]
  coeffs_out[1,7] <- fall.noJul.coeffs[2,1] 
  coeffs_out[2,1] <- sp.Jul.coeffs[1,1]
  coeffs_out[2,2] <- sp.Jul.coeffs[2,1]
  coeffs_out[2,4] <- sp.Jul.coeffs[3,1]
  coeffs_out[2,6] <- fall.Jul.coeffs[1,1]
  coeffs_out[2,7] <- fall.Jul.coeffs[2,1]                                
  coeffs_out[2,9] <- fall.Jul.coeffs[3,1]
  coeffs_out[3,1] <- sp.poly.coeffs[1,1]
  coeffs_out[3,2] <- sp.poly.coeffs[2,1]
  coeffs_out[3,3] <- sp.poly.coeffs[3,1]
  coeffs_out[3,6] <- fall.poly.coeffs[1,1]
  coeffs_out[3,7] <- fall.poly.coeffs[2,1]
  coeffs_out[3,8] <- fall.poly.coeffs[3,1]
  coeffs_out[4,1] <- sp.elev.coeffs[1,1]
  coeffs_out[4,2] <- sp.elev.coeffs[2,1]
  coeffs_out[4,4] <- sp.elev.coeffs[3,1]
  coeffs_out[4,5] <- sp.elev.coeffs[4,1]
  coeffs_out[4,6] <- fall.elev.coeffs[1,1]
  coeffs_out[4,7] <- fall.elev.coeffs[2,1]
  coeffs_out[4,9] <- fall.elev.coeffs[3,1]
  coeffs_out[4,10] <- fall.elev.coeffs[4,1]
  coeffs_out[5,1] <- sp.all.coeffs[1,1]
  coeffs_out[5,2] <- sp.all.coeffs[2,1]
  coeffs_out[5,3] <- sp.all.coeffs[3,1]
  coeffs_out[5,4] <- sp.all.coeffs[4,1]
  coeffs_out[5,5] <- sp.all.coeffs[5,1]
  coeffs_out[5,6] <- fall.all.coeffs[1,1]
  coeffs_out[5,7] <- fall.all.coeffs[2,1]
  coeffs_out[5,8] <- fall.all.coeffs[3,1]
  coeffs_out[5,9] <- fall.all.coeffs[4,1]
  coeffs_out[5,10] <- fall.all.coeffs[5,1]
  
  
  
  metrics_out <- round(metrics_out, digits = 2)
  coeffs_out <- round(coeffs_out, digits = 4)
  


write.table(x=coeffs_out, append=F,row.names=T, file = "2013_mod_coeffs_Mn.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = "2013_mod_metrics_Mn.csv", sep = ",", col.names=T)

write.table(x=NoNA.xyz, append=F,row.names=F, file = "2013_model_data_Mn.csv", sep = ",", col.names=T)


#####################################

#Below are some scripts that aren't used all the time, but may be useful for summary stuff




####################################################################################
#This part reads in the training dataset for a year and generates RMSEs for site-specific models
#####################################################################################

RcaID.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/RCA_DCE_unique.csv")
Log.in <- read.csv("JohnDay_8Day_2013.csv")
rca_elev <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/RCA_zone_min_elev.csv")

pathname <-"C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/"
RcaID.in <- merge(RcaID.in, rca_elev, by.x = "RASTERVALU", by.y = "RASTERVALU", all.x = TRUE, all.y = FALSE)

setwd("D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/31Jan13modelling/")
LST.in <- read.dbf("D:/Dropbox/work/research/CHaMP/CHaMP_data/LST/LST13_JD_RCA_smoothed.dbf")

merged.data <- merge(RcaID.in, Log.in, by.x = "DCE", by.y = "DceName", all.x = TRUE, all.y = TRUE)

HUCs.in <- read.csv("C:/Dropbox/work/research/Steelhead/spatial_data/coverages/RCA_HUCs.csv")

merged.data <- merge(merged.rca.hucs, Log.in, by.x = "SiteID", by.y = "SiteName", all.x = TRUE, all.y = FALSE)

ind <- apply(merged.data, 1, function(x) !any(is.na(x)))
Log.data <- merged.data[ind,]
Log.data$JulianDate <- as.numeric(Log.data$JulianDate)
Log.data <- orderBy(~RASTERVALU+JulianDate, Log.data)

colnames(Log.data)[6] <- "Mean"

merged_data_out <- data.frame (Juldate=NULL, DCE=NULL, RASTERVALU=NULL, Mean=NULL, y=NULL, LST=NULL)



colnames(LST.in)<-gsub("X", "", colnames(LST.in))
LST.in.names<-colnames(LST.in)
LST.in.names<- as.numeric(LST.in.names)
LST.in.names <- order(LST.in.names)
LST.in[,1:365] <- LST.in[,1:365] - 273.15


Rastervalu <- Log.data$OBJECTID
unique.IDs <- unique(Rastervalu)
merged_data_out <- data.frame (Juldate=NULL, DCE=NULL, RASTERVALU=NULL, Mean=NULL, y=NULL, LST=NULL)

coeffs_out <- data.frame(Site=numeric(length(unique.IDs)), sp_wJul_Int=numeric(length(unique.IDs)), sp_wJul_bLST=numeric(length(unique.IDs)), sp_wJul_bJul=numeric(length(unique.IDs)), fall_wJul_Int=numeric(length(unique.IDs)), fall_wJul_bLST=numeric(length(unique.IDs)), fall_wJul_bJul=numeric(length(unique.IDs)), sp_wJul_Int=numeric(length(unique.IDs)), sp_wJul_bLST=numeric(length(unique.IDs)), fall_wJul_Int=numeric(length(unique.IDs)), fall_wJul_bLST=numeric(length(unique.IDs)))                                                                                                                                                 

stats_out <- data.frame(Site=numeric(length(unique.IDs)), sp_wJul_rmse=numeric(length(unique.IDs)), sp_wJul_r2=numeric(length(unique.IDs)), N_spring=numeric(length(unique.IDs)), fall_wJul_rmse=numeric(length(unique.IDs)), fall_wJul_r2=numeric(length(unique.IDs)), sp_noJul_rmse=numeric(length(unique.IDs)), sp_noJul_r2=numeric(length(unique.IDs)), fall_noJul_rmse=numeric(length(unique.IDs)), fall_noJul_r2=numeric(length(unique.IDs)), N_fall=numeric(length(unique.IDs)))

for (i in 1:length(unique.IDs))
    {
    
    

     ind <- apply(mergeLogLST, 1, function(x) !any(is.na(x)))
     NoNA.xyz <- mergeLogLST[ind,]
     NoNA.xyz <- orderBy(~JulianDate, NoNA.xyz)

     maxrow <- which.max(NoNA.xyz$Mean)

     y.spring <- NoNA.xyz$Mean[1:maxrow]
     x.spring <- NoNA.xyz$LST[1:maxrow]
     z.spring <- NoNA.xyz$JulianDate[1:maxrow]
     dce.spring <- NoNA.xyz$RASTERVALU[1:maxrow]
  for (j in 1:6)
    {

     stream.log <- LST.Log.HUC[LST.Log.HUC$stream_ord == j,]
     
     
     ind <- apply(stream.log, 1, function(x) !any(is.na(x)))
     NoNA.xyz <- stream.log[ind,]
     NoNA.xyz <- orderBy(~JulianDate, NoNA.xyz)
     
              
     maxrow <- which.max(NoNA.xyz$Avg8Day)

     y.spring <- NoNA.xyz$Avg8Day[1:maxrow]
     x.spring <- NoNA.xyz$LST[1:maxrow]
     
     
     
     
     lst.sum.sp.huc <- summary(lm(y.spring ~ x.spring + I(x.spring^2)))
     pred.y.spring <- predict(lm(y.spring ~ x.spring + I(x.spring^2)))
     sp.poly.coeffs <- coefficients(lst.sum.sp.huc)
     sp.resids <- residuals(lst.sum.sp.huc)
     write.table (x=sp.resids,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/resids_sp_linear_Mn.csv",sep = ",", col.names=F)       
     
     
     pred.out <- matrix(nrow=length(y.spring), ncol=3)
     pred.out[,1] <- y.spring
     pred.out[,2] <- pred.y.spring 
     pred.out[,3] <- as.character(NoNA.xyz$SiteName[1:maxrow])
     write.table (x=pred.out,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/Pred_sp_poly_smooth.csv",sep = ",", col.names=F)  
     remove(pred.out)
     
     y.fall <- NoNA.xyz$Avg8Day[maxrow:nrow(NoNA.xyz)]
     x.fall <- NoNA.xyz$LST[maxrow:nrow(NoNA.xyz)]
     

   
     lst.sum.fall.huc<- summary(lm(y.fall ~ x.fall + I(x.fall^2)))
     pred.y.fall <- predict(lm(y.fall ~ x.fall + I(x.fall^2)))
     fall.poly.coeffs <- coefficients(lst.sum.fall.huc)
     fall.resids <- residuals(lst.sum.fall.huc)
     write.table (x=fall.resids,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/resids_fall_linear_Mn.csv",sep = ",", col.names=F)       
     
     pred.out <- matrix(nrow=length(y.fall), ncol=3)
     pred.out[,1] <- y.fall
     pred.out[,2] <- pred.y.fall  
     pred.out[,3] <- as.character(NoNA.xyz$SiteName[maxrow:nrow(NoNA.xyz)])
     write.table (x=pred.out,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/Pred_fall_poly_smooth.csv",sep = ",", col.names=F)  
     remove(pred.out)
            
     coeffs_out[j,1] <- j
     coeffs_out[j,2] <- sp.poly.coeffs[1,1]
     coeffs_out[j,3] <- sp.poly.coeffs[2,1]
     coeffs_out[j,4] <- sp.poly.coeffs[3,1]
     coeffs_out[j,5] <- fall.poly.coeffs[1,1]
     coeffs_out[j,6] <- fall.poly.coeffs[2,1]
     coeffs_out[j,7] <- fall.poly.coeffs[3,1]
     
     
     stats_out[j,1] <- j
     stats_out[j,2] <- lst.sum.sp.huc$sigma
     stats_out[j,3] <- lst.sum.sp.huc$adj.r.squared
     stats_out[j,4] <- lst.sum.fall.huc$sigma
     stats_out[j,5] <- lst.sum.fall.huc$adj.r.squared
     

     
    }


write.table (x=coeffs_out,append=F,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/2009_indiv_model_coeffs.csv",sep = ",", col.names=T)
 
write.table (x=stats_out,append=F,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/2009_indiv_model_stats.csv",sep = ",", col.names=T)
write.table (x=stats_out,append=F,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/metrics_output/2009_indiv_model_stats.csv",sep = ",", col.names=T)
 
rm(list = ls())

####################################################################################
#This part reads in the validation and/or training dataset for a year and compares it site-by-site to the predicted values
#####################################################################################


RcaID.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/RCA_DCE_unique.csv")
Log.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/training_sets/All_sites_2009_logger_train.csv")
rca_elev <- read.csv("C:/Dropbox/work/research/Steelhead/spatial_data/RCA_zone_min_elev.csv")

RcaID.in <- merge(RcaID.in, rca_elev, by.x = "RASTERVALU", by.y = "RASTERVALU", all.x = TRUE, all.y = FALSE)

merged.data <- merge(RcaID.in, Log.in, by.x = "SiteName", by.y = "SiteName", all.x = TRUE, all.y = TRUE)
ind <- apply(merged.data, 1, function(x) !any(is.na(x)))
Log.data <- merged.data[ind,]
Log.data$JulianDate <- as.numeric(Log.data$JulianDate)
Log.data <- orderBy(~GRIDCODE+JulianDate, Log.data)

#colnames(Log.data)[6] <- "Mean"

#use appropriate read statement
predT.in <- read.csv("D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/Pred_fall_poly_smooth.csv")
#predT.in <- tLogPred.out        



GRIDCODE <- Log.data$GRIDCODE
unique.IDs <- unique(GRIDCODE)
coeffs_out <- data.frame(GRIDCODE=(length(unique.IDs)), r2=numeric(length(unique.IDs)), RMSE=numeric(length(unique.IDs)))
mergeLogLST <- data.frame (mup = NULL)

for (i in 1:length(unique.IDs))
    {
    Log.site <- Log.data[Log.data$GRIDCODE == unique.IDs[i],]
    
    predT.day <- matrix(nrow=46, ncol=2)
    predT.day[,1] <- as.numeric(colnames(predT.in)[1:46])
    predT.day[,2] <- unlist(predT.in[predT.in$GRID_CODE == unique.IDs[i],1:46])
    predT.day <- data.frame(predT.day)
    colnames(predT.day) <- c("JulDay", "predT")
    mergeLogLST <- merge(Log.site, predT.day, by.x = "JulianDate", by.y = "JulDay", all.x = TRUE, all.y = TRUE)
    
    ind <- apply(mergeLogLST, 1, function(x) !any(is.na(x)))
    NoNA.xyz <- mergeLogLST[ind,]
    
       
    y <- NoNA.xyz$Avg8Day
    x <- NoNA.xyz$predT
    
    

    lst.sum <- summary(lm(y ~ x))
    pred.y <- predict(lm(y ~ x ))
    coeff <- coefficients(lst.sum)
       
    pred.out <- matrix(nrow=length(y), ncol=3)
    pred.out[,1] <- y
    pred.out[,2] <- x 
    pred.out[,3] <- y - x 
    
    colnames(pred.out) <- c("Logger", "PredT", "Resids")
     
    write.table (x=pred.out,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Asotin/18Aug14modeling/Y_predY_2013_noJul_train.csv",sep = ",", col.names=F)
    write.table (x=mergeLogLST,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Asotin/18Aug14modeling/2013_noJul_train_preds.csv",sep = ",", col.names=F)
    
    rmse <- (y - x)^2   
    coeffs_out[i,1] <- unique.IDs[i]
    coeffs_out[i,2] <- lst.sum$adj.r.squared
    coeffs_out[i,3] <- sqrt(mean(rmse))
        
       
    remove(lst.sum)
    remove(pred.y)
    remove(coeff)
    remove(pred.out)
    remove(x)
    remove(y)
    
    }
    
coeffs_out$r2 <- round(coeffs_out$r2, digits = 2)
coeffs_out$RMSE <- round(coeffs_out$RMSE, digits = 2)        

write.table (x=coeffs_out,append=F,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/RMSE_2009_noJul_train_site_by_predicted.csv",sep = ",", col.names=T)


rm(list = ls())

####################################################################################
#This part reads in the validation and/or training dataset for a year and compares those data globally to the predicted values
#####################################################################################


RcaID.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/RCA_DCE_unique.csv")
Log.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/training_sets/All_sites_2009_logger_train.csv")
rca_elev <- read.csv("C:/Dropbox/work/research/Steelhead/spatial_data/RCA_zone_min_elev.csv")

RcaID.in <- merge(RcaID.in, rca_elev, by.x = "RASTERVALU", by.y = "RASTERVALU", all.x = TRUE, all.y = FALSE)

merged.data <- merge(RcaID.in, Log.in, by.x = "DCE", by.y = "DceName", all.x = TRUE, all.y = TRUE)
ind <- apply(merged.data, 1, function(x) !any(is.na(x)))
Log.data <- merged.data[ind,]
Log.data$JulianDate <- as.numeric(Log.data$JulianDate)
Log.data <- orderBy(~RASTERVALU+JulianDate, Log.data)

colnames(Log.data)[6] <- "Mean"

#use appropriate read statement
predT.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/Predicted_RCAs/predt2001_noJul.csv")
#predT.in <- tLogPred.out        



Rastervalu <- Log.data$RASTERVALU
unique.IDs <- unique(Rastervalu)
coeffs_out <- data.frame (Data_set=(1), r2=(1), RMSE=(1))
merged_data_out <- data.frame (Juldate=NULL, DCE=NULL, RASTERVALU=NULL, Elev=NULL, SampleDate=NULL, Year=NULL, Mean=NULL, predT=NULL)  

for (i in 1:length(unique.IDs))
    {
    Log.site <- Log.data[Log.data$OBJECTID == unique.IDs[i],]
    
    predT.day <- matrix(nrow=length(Log.site), ncol=3)
    predT.day[,1] <- as.numeric(colnames(predT.in)[2:47]
    predT.day[,2] <- unlist(predT.in[predT.in$OBJECTID == unique.IDs[i],2:47])
    predT.day <- data.frame(predT.day)
    colnames(predT.day) <- c("JulDay", "predT")
    mergeLogLST <- merge(Log.site, predT.day, by.x = "JulianDate", by.y = "JulDay", all.x = TRUE, all.y = TRUE)
    
    ind <- apply(mergeLogLST, 1, function(x) !any(is.na(x)))
    NoNA.xyz <- mergeLogLST[ind,]
    
    merged_data_out <- rbind(merged_data_out, NoNA.xyz)
    }
 

y <- merged_data_out$Mean
x <- merged_data_out$predT

lst.sum <- summary(lm(y ~ x))
pred.y <- predict(lm(y ~ x ))
coeff <- coefficients(lst.sum)
   
pred.out <- matrix(nrow=length(y), ncol=4)
pred.out[,1] <- y
pred.out[,2] <- x 
pred.out[,3] <- y - x 
pred.out[,4] <- "2009"
colnames(pred.out) <- c("Logger", "PredT", "Resids", "Year")
 
write.table (x=pred.out,append=F,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/Y_predY_global_2009_noJul_train.csv",sep = ",", col.names=T)

write.table (x=merged_data_out,append=T,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/2009_noJul_train_merged_data.csv",sep = ",", col.names=F)

rmse <- (y - x)^2   
coeffs_out$Data_set <- "2009_noJul_train"
coeffs_out$r2 <- round(lst.sum$adj.r.squared, digits = 2)
coeffs_out$RMSE <- round(sqrt(mean(rmse)), digits = 2)
        
       
write.table (x=coeffs_out,append=F,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/RMSE_global_data_on_predicted_2009.csv",sep = ",", col.names=F)

write.table (x=coeffs_out,append=T,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/metrics_output/RMSE_global_data_on_predicted.csv",sep = ",", col.names=F)

     
write.table (x=merged.data,append=F,row.names=F,file="C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/global_all_data_2009_noJul_train.csv",sep = ",", col.names=T)   


