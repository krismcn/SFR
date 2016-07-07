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

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/shapes/")

  LST.in <- read.dbf("USal_13_LST.dbf")
  PointID <- LST.in[,1]
  LST.in <- LST.in[,3:48]
  LST.in[LST.in<0.1] = NA
  
  
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
  LST.out$PointID <- PointID
  

  write.dbf(LST.out, file = "LST13_USal_interp.dbf")

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################

  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/")
  out_Preds <- read.dbf("LST13_USal_interp.dbf") 
  
  weights <- read.dbf("USal_rca_area_wgts.dbf")
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
          daily <- out_Preds[out_Preds$PointID %in% pixels, j]
          zonal_mn <- weighted.mean(daily, wgts, na.rm = TRUE)
          rca_zonal[l,j] <- zonal_mn
          }
      l <- l+1
      }
  

colnames(rca_zonal)[1:46] <- colnames(out_Preds)[1:46]
rca_zonal <- as.data.frame(rca_zonal)
rca_zonal$RCAID <- rcas


write.dbf(rca_zonal, file = "LST13_USal_RCA.dbf")



#################################################################
#Logger summary
#Fills in date gaps with NAs and generates 8-day means 
#################################################################
  
  ID.in <- read.csv("USal_site_elev.csv")
  LST.in <- read.dbf("LST13_USal_RCA.dbf")

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/USalmon/2013/")

##############
# NorWeST data
##############
  Log.in <- read.csv("USal_2013_NorWest_logger_data.csv")
  
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    #specifies those dates to be numeric in order to drop the leading zeros
  #newnamesnum <- as.numeric(colnames(LST.in)[1:46]);
  colnames(LST.in)[1:46] <- newnamesnum;
  
  dates <- matrix(nrow=365, ncol=1)
  dates[,1] <- 1:365
  colnames(dates) <- c("Days")
  
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)
  Log.8Day.out <- data.frame (mup = NULL)
  
    for (i in SiteID) 
      { 
        Log.site <- Log.in[Log.in$SiteName == i,]
        full.year <- merge(dates, Log.site, by.x = "Days", by.y = "JulDay", all.x = TRUE)
        eightday <- rollmean(full.year$DailyMaximum, 8, fill = NA, align = "left")
        eightday <- as.matrix(eightday)
        full.year$Mx8D <- eightday
        full.year$Mx8D[361] <- mean(full.year$DailyMaximum[361:365])
        Log.8Day.out <- rbind(Log.8Day.out, full.year)
        
      }

  write.csv(x = Log.8Day.out, file = "NorWest_2013_8Day_logger_data.csv", row.names = F)

##############
# CHaMP data
##############

  Log.in <- read.csv("USal_2013_CHaMP_logger_data.csv")
  
  SiteID <- unique(Log.in$SiteName)
  SiteID <- as.matrix(SiteID)
  Log.8Day.out <- Log.in


##############
# All data
##############

Log.in <- read.csv("All_data_2013_8D_logger.csv")

SiteID <- unique(Log.in$SiteName)
SiteID <- as.matrix(SiteID)
Log.8Day.out <- Log.in

#################################################################
#Logger prediction modeling part
#Generates models for both LST-only and LST + Julian Day datasets
#################################################################
  
  LST.Log.out <- data.frame (mup = NULL)

  for (i in SiteID) 
    { 
      Log.site <- Log.8Day.out[Log.8Day.out$SiteName  == i,]
      Log.site <- as.data.frame(Log.site)
      
      
      RCAID <- ID.in$rca_id[ID.in$SiteName == i]
      Elev <- ID.in$ELEV[ID.in$SiteName == i]
      
      LST.site <- matrix(ncol=3, nrow=46)
      LST.site[,1] <- as.numeric(unlist(colnames(LST.in)[1:46]))
      LST.site[,2] <- unlist(LST.in[RCAID,1:46])
      LST.site <- data.frame(LST.site)
      colnames(LST.site) <- c("JulDay", "LST", "Elev")
      LST.site[3] <- Elev
      LST.Log.site <- merge(LST.site, Log.site, by.x = "JulDay", by.y = "JulDay", all.x=TRUE, all.y = FALSE)
      
      LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
    }
    
coeffs_out <- data.frame(sp_Int=numeric(5), sp_bLST=numeric(5), sp_bLST2=numeric(5), sp_bJul=numeric(1), sp_bElev=numeric(5))
metrics_out <- data.frame(PRESS_sp=numeric(5), p2_sp=numeric(5), r2_sp=numeric(5), RMSEP_sp=numeric(5), RMSE_sp=numeric(5), AIC_sp=numeric(5))
rownames(metrics_out) <- c("noJul", "Jul", "Poly", "Elev", "All")
rownames(coeffs_out) <- c("noJul", "Jul", "Poly", "Elev", "All")

    ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
    NoNA.xyz <- LST.Log.out[ind,]

write.csv(x=NoNA.xyz, file="USal_2013_8Day_CHaMP_model_data.csv")
  
    NoNA.xyz <- NoNA.xyz[,c(7, 2, 1, 3, 4)]
    colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
    NoNA.xyz$SiteName <- as.character(NoNA.xyz$SiteName)
    NoNA.xyz <- orderBy(~z, NoNA.xyz)

####spring#######

  maxrow <- which.max(NoNA.xyz$y)
  y <- NoNA.xyz$y[1:maxrow]
  x <- NoNA.xyz$x[1:maxrow]
  z <- NoNA.xyz$z[1:maxrow]
  e <- NoNA.xyz$e[1:maxrow]
  plot(x,y)

####fall#######

  maxrow <- which.max(NoNA.xyz$y)
  y <- NoNA.xyz$y[maxrow:nrow(NoNA.xyz)]
  x <- NoNA.xyz$x[maxrow:nrow(NoNA.xyz)]
  z <- NoNA.xyz$z[maxrow:nrow(NoNA.xyz)]
  e <- NoNA.xyz$e[maxrow:nrow(NoNA.xyz)]
  plot(x,y)

########################
    summer <- subset(NoNA.xyz, z > 181 & z < 258)
    SiteID <- as.matrix(unique(summer$SiteName))
    SiteID <- as.matrix(as.character(SiteID))
    ID.table <- as.data.frame(table(summer$SiteName))  
    SiteID.sum <- as.matrix(ID.table[ID.table$Freq > 4, "Var1"])
    
    y <- summer$y
    x <- summer$x
    z <- summer$z
    e <- summer$e
    plot(x, y)

###########################    
    y <- NoNA.xyz$y
    x <- NoNA.xyz$x
    z <- NoNA.xyz$z
    e <- NoNA.xyz$e
    plot(x, y)

############################
# single model
###########################

    noJul_mod <- lm(y ~ x + I(x^2) + e)
    noJul_summ_sp <- summary(lm(y ~ x + I(x^2) + e))
    noJul_plsr_mod_sp <- plsr(y ~ x + I(x^2) + e, validation = "LOO")
    noJul_RMSEP_sp <- RMSEP(noJul_plsr_mod_sp)
    sp.noJul.coeffs <- coefficients(noJul_summ_sp)
    
    Jul_mod <- lm(y ~ x + z)
    Jul_summ_sp <- summary(lm(y ~ x + z))
    Jul_plsr_mod_sp <- plsr(y ~ x + z, validation = "LOO")
    Jul_RMSEP_sp <- RMSEP(Jul_plsr_mod_sp)
    sp.Jul.coeffs <- coefficients(Jul_summ_sp)
    
    poly_mod <- lm(y ~ x + I(x^2))
    poly_summ_sp <- summary(lm(y ~ x + I(x^2) +z))
    poly_plsr_mod_sp <- plsr(y ~ x + I(x^2) + z , validation = "LOO")
    poly_RMSEP_sp <- RMSEP(poly_plsr_mod_sp)
    sp.poly.coeffs <- coefficients(poly_summ_sp)
    
    elev_mod_sp <- lm(y ~ x + z + e)
    elev_summ_sp <- summary(lm(y ~ x + z + e))     
    elev_plsr_mod_sp <- plsr(y ~ x + z + e, validation = "LOO")
    elev_RMSEP_sp <- RMSEP(elev_plsr_mod_sp)
    sp.elev.coeffs <- coefficients(elev_summ_sp)
    
    all_mod_sp <- lm(y ~ x + I(x^2) + z + e)
    all_summ_sp <- summary(lm(y ~ x + I(x^2) + z + e))   
    all_plsr_mod_sp <- plsr(y ~ x + I(x^2) + z + e, validation = "LOO")
    all_RMSEP_sp <- RMSEP(all_plsr_mod_sp)    
    sp.all.coeffs <- coefficients(all_summ_sp)
    
    pressstat_Jul <- PRESS(Jul_mod, verbose = "FALSE")
    pressstat_noJul <- PRESS(noJul_mod, verbose = "FALSE")
    pressstat_poly <- PRESS(poly_mod, verbose = "FALSE")
    pressstat_elev <- PRESS(elev_mod_sp, verbose = "FALSE")
    pressstat_all <- PRESS(all_mod_sp, verbose = "FALSE")
    
    pressstat_Jul$stat
    pressstat_noJul$stat
    pressstat_poly$stat
    pressstat_elev$stat
    pressstat_all$stat
    
    AIC(noJul_mod)
    AIC(Jul_mod)
    AIC(poly_mod)
    AIC(elev_mod_sp)
    AIC(all_mod_sp)
    
    # output ------------------------------------------------------------------
    
    
    metrics_out[1,1] <- pressstat_noJul$stat
    metrics_out[1,2] <- pressstat_noJul$P.square
    metrics_out[1,3] <- noJul_summ_sp$adj.r.squared     
    metrics_out[1,4] <- noJul_RMSEP_sp[[1]][[4]]
    metrics_out[1,5] <- noJul_summ_sp$sigma
    metrics_out[1,6] <- AIC(noJul_mod)
    metrics_out[2,1] <- pressstat_Jul$stat
    metrics_out[2,2] <- pressstat_Jul$P.square
    metrics_out[2,3] <- Jul_summ_sp$adj.r.squared     
    metrics_out[2,4] <- Jul_RMSEP_sp[[1]][[4]]
    metrics_out[2,5] <- Jul_summ_sp$sigma
    metrics_out[2,6] <- AIC(Jul_mod)
    metrics_out[3,1] <- pressstat_poly$stat
    metrics_out[3,2] <- pressstat_poly$P.square
    metrics_out[3,3] <- poly_summ_sp$adj.r.squared     
    metrics_out[3,4] <- poly_RMSEP_sp[[1]][[4]]
    metrics_out[3,5] <- poly_summ_sp$sigma
    metrics_out[3,6] <- AIC(poly_mod)
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

# output ------------------------------------------------------------------


    if (metrics_out[1,4] < metrics_out[2,4] & metrics_out[1,4] < metrics_out[3,4] & metrics_out[1,4] < metrics_out[4,4] & metrics_out[1,4] < metrics_out[5,4])
    {pred.y<- predict(noJul_mod); sp_mod <- noJul_mod}
    metrics_out[1,4] < metrics_out[2,4] & metrics_out[1,4] < metrics_out[3,4] & metrics_out[1,4] < metrics_out[4,4] & metrics_out[1,4] < metrics_out[5,4]
    
    if (metrics_out[2,4] < metrics_out[1,4] & metrics_out[2,4] < metrics_out[3,4] & metrics_out[2,4] < metrics_out[4,4] & metrics_out[2,4] < metrics_out[5,4])
    {pred.y.spring <- predict(Jul_mod); sp_mod <- Jul_mod}
    metrics_out[2,4] < metrics_out[1,4] & metrics_out[2,4] < metrics_out[3,4] & metrics_out[2,4] < metrics_out[4,4] & metrics_out[2,4] < metrics_out[5,4]      
    
    if (metrics_out[3,4] < metrics_out[2,4] & metrics_out[3,4] < metrics_out[1,4] & metrics_out[3,4] < metrics_out[4,4] & metrics_out[3,4] < metrics_out[5,4])
    {pred.y.spring <- predict(poly_mod); sp_mod <- poly_mod}
    metrics_out[3,4] < metrics_out[2,4] & metrics_out[3,4] < metrics_out[1,4] & metrics_out[3,4] < metrics_out[4,4] & metrics_out[3,4] < metrics_out[5,4]
    
    if (metrics_out[4,4] < metrics_out[2,4] & metrics_out[4,4] < metrics_out[1,4] & metrics_out[4,4] < metrics_out[3,4] & metrics_out[4,4] < metrics_out[5,4])
    {pred.y.spring <- predict(elev_mod_sp); sp_mod <- elev_mod_sp}
    metrics_out[4,4] < metrics_out[2,4] & metrics_out[4,4] < metrics_out[1,4] & metrics_out[4,4] < metrics_out[3,4] & metrics_out[4,4] < metrics_out[5,4]      
    
    if (metrics_out[5,4] < metrics_out[2,4] & metrics_out[5,4] < metrics_out[1,4] & metrics_out[5,4] < metrics_out[3,4] & metrics_out[5,4] < metrics_out[4,4])
    {pred.y.spring <- predict(all_mod_sp); sp_mod <- all_mod_sp}
    metrics_out[5,4] < metrics_out[2,4] & metrics_out[5,4] < metrics_out[1,4] & metrics_out[5,4] < metrics_out[3,4] & metrics_out[5,4] < metrics_out[4,4]      
    
    pred.y.spring[pred.y.spring<0.5] = 0
    pred.out <- matrix(nrow=length(y), ncol=4)
    colnames(pred.out) <- c("Y", "PredY", "SiteName", "Year")
    pred.out[,1] <- y
    pred.out[,2] <- pred.y.spring
    pred.out[,3] <- as.character(NoNA.xyz$SiteName)
    pred.out[,4] <- "USal"
    plot(pred.out[,1], pred.out[,2])
    
    
    write.table (x=pred.out,append=F,row.names=F,file="2013_pred_v_y_USal.csv",sep = ",", col.names=T)  
    
    
    
    coeffs_out[1,1] <- sp.noJul.coeffs[1,1]
    coeffs_out[1,2] <- sp.noJul.coeffs[2,1]
    coeffs_out[1,3] <- sp.noJul.coeffs[3,1]
    coeffs_out[1,5] <- sp.noJul.coeffs[4,1]
    
    coeffs_out[2,1] <- sp.Jul.coeffs[1,1]
    coeffs_out[2,2] <- sp.Jul.coeffs[2,1]
    coeffs_out[2,4] <- sp.Jul.coeffs[3,1]
    
    coeffs_out[3,1] <- sp.poly.coeffs[1,1]
    coeffs_out[3,2] <- sp.poly.coeffs[2,1]
    coeffs_out[3,3] <- sp.poly.coeffs[3,1]
    
    coeffs_out[4,1] <- sp.elev.coeffs[1,1]
    coeffs_out[4,2] <- sp.elev.coeffs[2,1]
    coeffs_out[4,4] <- sp.elev.coeffs[3,1]
    coeffs_out[4,5] <- sp.elev.coeffs[4,1]
    
    coeffs_out[5,1] <- sp.all.coeffs[1,1]
    coeffs_out[5,2] <- sp.all.coeffs[2,1]
    coeffs_out[5,3] <- sp.all.coeffs[3,1]
    coeffs_out[5,4] <- sp.all.coeffs[4,1]
    coeffs_out[5,5] <- sp.all.coeffs[5,1]
    
    
    
    
    metrics_out <- round(metrics_out, digits = 2)
    coeffs_out <- round(coeffs_out, digits = 4)
    
    rownames(metrics_out) <- c("noJul", "Jul", "Poly", "Elev", "All")
    rownames(coeffs_out) <- c("noJul", "Jul", "Poly", "Elev", "All")
    
    
    write.table(x=coeffs_out, append=F,row.names=F, file = "All_data_2013_mod_coeffs_Mn.csv", sep = ",", col.names=T)
    
    write.table(x=metrics_out, append=F,row.names=F, file = "All_data_2013_mod_metrics_Mn.csv", sep = ",", col.names=T)

#####################################
# spring/fall
######################################

    setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/2013/")

    coeffs_out <- data.frame(Int=numeric(2), bLST=numeric(2), bLST2=numeric(2), bJul=numeric(2), bElev=numeric(2))
    metrics_out <- data.frame(PRESS=numeric(2), p2=numeric(2), r2=numeric(2), RMSEP=numeric(2), RMSE=numeric(2))
    rownames(metrics_out) <- c("Spring", "Fall")
    rownames(coeffs_out) <- c("Spring", "Fall")

    maxrow <- which.max(NoNA.xyz$y)
    data.sp <- NoNA.xyz[1:maxrow,]
    data.fall <- NoNA.xyz[maxrow:nrow(NoNA.xyz),]
    
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

    pred.out <- matrix(nrow=length(pred.y), ncol=5)
    pred.out[,1] <- y
    pred.out[,2] <- pred.y
    pred.out[,3] <- z
    pred.out[,4] <- "spring"
    pred.out[,5] <- "2013"
    colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year")

    RSS <- sum(pressstat_sum$residuals^2)
    TSS <- sum((y - mean(y))^2)
    MSS <- TSS - RSS
    p2 <- MSS/TSS
    
   
    write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_USal_2013_yr.csv",sep = ",", col.names=F)  
    
    
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


#__________________________

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
    
    pred.out <- matrix(nrow=length(pred.y), ncol=5)
    pred.out[,1] <- y
    pred.out[,2] <- pred.y
    pred.out[,3] <- z
    pred.out[,4] <- "fall"
    pred.out[,5] <- "2013"
    
    write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_USal_2013_yr.csv",sep = ",", col.names=F)  

    RSS <- sum(pressstat_sum$residuals^2)
    TSS <- sum((y - mean(y))^2)
    MSS <- TSS - RSS
    p2 <- MSS/TSS
    
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


  metrics_out <- round(metrics_out, digits = 2)
  coeffs_out <- round(coeffs_out, digits = 4)

write.table(x=coeffs_out, append=F,row.names=T, file = "All_data_2013_mod_coeffs_Mx.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=T, file = "All_data_2013_mod_metrics_Mx.csv", sep = ",", col.names=T)

  #______________________
  # went back and calc-ed the model metrics for sping & fall together:
    
    setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/2013/")

    Log.in <- read.csv("USal_2013_8Day_CHaMP_model_data.csv")
    NoNA.xyz <- Log.in[,c(8, 3, 2, 4, 5)]
    colnames(NoNA.xyz) <- c("y", "x", "z", "e", "SiteName")
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


########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################


  
  LST.in <- read.dbf("D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/USalmon/LST13_USal_RCA.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  namesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5))
  colnames(LST.in)[1:46] <- namesnum;
  
  elev.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/USalmon/USal_rca_elev.csv")
  LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "rca_id")
  
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
  LogPred.out$Basin_RCA <- paste0("USal_", LogPred.out$RCAID)
  namesnum <- as.numeric(colnames(LogPred.out[1:12]))
  
  names.out <- sprintf("Tmx_12_%03d", namesnum)
  colnames(LogPred.out)[1:12] <- names.out[1:12]
  
  write.dbf(LogPred.out, file = "predt2013_USal_8D_Max_summer.dbf") 

#________________________________________________________
# summary for the summer max EPs
# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/EP_temp/USalmon/"
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

  Au.in <- read.csv("USal_RCA_AUs_IP.csv")
  ST_rcas <- Au.in[Au.in$STHD_IP == "Steelhead",]
  CH_rcas <- Au.in[Au.in$CHNK_IP == "Chinook",]
  
  setwd(paste0(mainPath,yearPath))
  
  Max.in <- read.dbf("predt2013_USal_8D_Max_summer.dbf") 
  
  Max.au <- merge(Max.in, Au.in, by.x = "RCAID", by.y = "rca_id")
  Max.au$AU_Code <- as.character(Max.au$AU_ID)
  Aus <- unique(Max.au$AU_Code)
  Aus <- Aus[1:2]
  
  Max.daily <- Max.in[, c(rep(3:8, each = 8), 13:15)]
  colnames(Max.daily)[1:48] <- 201:248
  
  #rcas <- unique(ST_rcas$rca_id)
  #rcas <- unique(CH_rcas$rca_id)
  rcas <- unique(Max.in$RCAID)
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

colnames(SumSumm.out) <- c("RCAID", "Pct12_2013","Pct13_2013", "Pct16_2013", "Pct18_2013", "Pct20_2013", "Pct22_2013",  "MxMx_2013", "sdMn_2013", "MnMx_2013")


write.dbf(SumSumm.out, file = "USal_2013_21Jul_31Aug_max_summary_All.dbf")
write.csv(SumSumm.out, file = "USal_2013_21Jul_31Aug_max_summary_All.csv", row.names = F)


# summary by AU -----------------------------------------------------------



  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "rca_id")
  colnames(SumSummAU) <- c("RCAID", "PctDays12","PctDays13", "PctDays16", "PctDays18", "PctDays20", "PctDays22",  "MxMx", "sdMn", "MnMx", "STHD_IP", "CHNK_IP", "AU_Code")
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
      SumAU <- SumSummAU[SumSummAU$AU_Code == Aus[i],]
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
write.csv(SumAU.out, file = "USal_AU_2013_21Jul_31Aug_ExPcnt_summary_All.csv", row.names = F)
