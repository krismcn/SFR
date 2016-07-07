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
setwd("D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat")

library(timeSeries)
library(lattice)
library(foreign)
#library(zoo)
library(doBy)
library(qpcR)

LST.in <- read.csv("Ent_LST_2012.csv", header = TRUE)
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

out_Preds <- read.dbf("LST12_Ent_interp.dbf") #use the appropriate read statement

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

write.dbf(rca_zonal, file = "LST12_Ent_RCA.dbf")



#################################################################
#Logger prediction modeling part
#Generates models for both LST-only and LST + Julian Day datasets
#################################################################
setwd("D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat")

ID.in <- read.csv("Ent_sites_elev.csv")
Log.in <- read.csv("Ent_temp_data_2012.csv")

#Log.in[Log.in<0] = 0.0
 
LST.in <- read.dbf("D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/LST12_Ent_RCA.dbf")
colnames(LST.in)<-gsub("X", "", colnames(LST.in))
#newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    #specifies those dates to be numeric in order to drop the leading zeros
newnamesnum <- as.numeric(colnames(LST.in)[1:46]);
colnames(LST.in)[1:46] <- newnamesnum;

#newnames <- substring (colnames(LST.in)[1:46],7,11); #clips the MODIS grid names down to the julian date
#newnamesnum <- as.numeric(newnames);    #specifies those dates to be numeric in order to drop the leading zeros
#colnames(LST.in)[1:46] <- newnamesnum;

#LST.zoo <- zoo(t(LST.in[,1:46]))
#LST.zoo.sm <- rollapply(LST.zoo, 2, (mean), by.column = TRUE, na.rm = T, align="right", partial = 1, fill = NA)
#LST.zoo.out <- as.data.frame(t(LST.zoo.sm))
#LST.zoo.out <- round(LST.zoo.out, digits = 2)
#colnames(LST.zoo.out) <- colnames(LST.in[1:46])
#rownames(LST.zoo.out) <- NULL

#LST.zoo.out$GRID_CODE <- GrPolID
#write.dbf(LST.zoo.out, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/LST12_Ent_kmpt_smoothed.dbf")

coeffs_out <- data.frame(sp_Jul_Int=numeric(4), sp_Jul_bLST=numeric(4), sp_Jul_bJul=numeric(4), fall_Jul_Int=numeric(4), fall_Jul_bLST=numeric(4), fall_Jul_bJul=numeric(4), sp_poly_Int=numeric(4), sp_poly_bLST=numeric(4), fall_poly_Int=numeric(4), fall_poly_bLST=numeric(4), sp_poly_bLST2=numeric(4), fall_poly_bLST2=numeric(4))
stats_out <- data.frame(sp_Jul_RMSE=numeric(4), sp_Jul_r2=numeric(4), fall_Jul_RMSE=numeric(4), fall_Jul_r2=numeric(4), sp_poly_RMSE=numeric(4), sp_poly_r2=numeric(4), fall_poly_RMSE=numeric(4), fall_poly_r2=numeric(4))
metrics_out <- data.frame(sp_Jul_PRESS=numeric(4), sp_Jul_p2=numeric(4), sp_noJul_PRESS = numeric(4), sp_noJul_p2=numeric(4), sp_poly_PRESS=numeric(4), sp_poly_p2=numeric(4), fall_Jul_PRESS=numeric(4), fall_Jul_p2=numeric(4), fall_noJul_PRESS = numeric(4), fall_noJul_p2=numeric(4), fall_poly_PRESS=numeric(4), fall_poly_p2=numeric(4))
Log.GRIDCODE <- merge(Log.in, ID.in, by.x = "SiteName", by.y = "SiteName", all.x = T, all.y = F)
LST.Log.out <- data.frame (mup = NULL)
SiteID <- unique(Log.in$SiteName)
SiteID <- as.matrix(SiteID)

for (i in SiteID) 
        { 
        Log.site <- Log.in[Log.in$SiteName == i,]
        Log.site <- as.data.frame(Log.site)
        
        
        GridID <- ID.in$Entiat_RCA[ID.in$SiteName == i]
        
        LST.site <- matrix(ncol=2, nrow=46)
        LST.site[,1] <- as.numeric(unlist(colnames(LST.in)[1:46]))
        LST.site[,2] <- unlist(LST.in[GridID,1:46])
        LST.site <- data.frame(LST.site)
        colnames(LST.site) <- c("JulDay", "LST")
        
        LST.Log.site <- merge(LST.site, Log.site, by.x = "JulDay", by.y = "JulDay", all.x=TRUE, all.y = FALSE)
        
        LST.Log.out <- rbind(LST.Log.out, LST.Log.site)
        }





     
     
     ind <- apply(LST.Log.out, 1, function(x) !any(is.na(x)))
     NoNA.xyz <- LST.Log.out[ind,]
     NoNA.xyz <- orderBy(~JulDay, NoNA.xyz)
     
              
     maxrow <- which.max(NoNA.xyz$Avg8Day)

     y.spring <- NoNA.xyz$Avg8Day[1:maxrow]
     x.spring <- NoNA.xyz$LST[1:maxrow]
     z.spring <- NoNA.xyz$JulDay[1:maxrow]
     
     noJul_mod <- lm(y.spring ~ x.spring)
     noJul_summ_sp <- summary(lm(y.spring ~ x.spring))


     Jul_mod <- lm(y.spring ~ x.spring + z.spring)
     Jul_summ_sp <- summary(lm(y.spring ~ x.spring + z.spring))
     poly_mod <- lm(y.spring ~ x.spring + I(x.spring^2))
     poly_summ_sp <- summary(lm(y.spring ~ x.spring + I(x.spring^2)))
     pred.y.spring.poly <- predict(lm(y.spring ~ x.spring + I(x.spring^2)))
     sp.poly.coeffs <- coefficients(poly_summ_sp)
     sp.Jul.coeffs <- coefficients(Jul_summ_sp)
     sp.resids <- residuals(poly_summ_sp)
     pressstat_Jul <- PRESS(Jul_mod)
     pressstat_Jul$stat
     pressstat_noJul <- PRESS(noJul_mod)
     pressstat_noJul$stat
     pressstat_poly <- PRESS(poly_mod)
     pressstat_poly$stat
     sp.resids <- residuals(poly_summ_sp)
     write.table (x=sp.resids,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/resids_sp_linear_Mn.csv",sep = ",", col.names=F)       
     
    metrics_out[1,1] <- pressstat_Jul$stat
    metrics_out[1,2] <- pressstat_Jul$P.square
    metrics_out[1,3] <- pressstat_noJul$stat
    metrics_out[1,4] <- pressstat_noJul$P.square
    metrics_out[1,5] <- pressstat_poly$stat
    metrics_out[1,6] <- pressstat_poly$P.square
     
     pred.out <- matrix(nrow=length(y.spring), ncol=3)
     pred.out[,1] <- y.spring
     pred.out[,2] <- pred.y.spring.poly 
     pred.out[,3] <- as.character(NoNA.xyz$SiteName[1:maxrow])
     write.table (x=pred.out,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/Pred_sp_poly_Mn.csv",sep = ",", col.names=F)  
     remove(pred.out)
     
     y.fall <- NoNA.xyz$Avg8Day[maxrow:nrow(NoNA.xyz)]
     x.fall <- NoNA.xyz$LST[maxrow:nrow(NoNA.xyz)]
     z.fall <- NoNA.xyz$JulDay[maxrow:nrow(NoNA.xyz)]

    noJul_mod <- lm(y.fall ~ x.fall)
    noJul_summ_sp <- summary(lm(y.fall ~ x.fall))

    Jul_mod <- lm(y.fall ~ x.fall + z.fall) 
    Jul_summ_fall<- summary(lm(y.fall ~ x.fall + z.fall))
     poly_summ_fall<- summary(lm(y.fall ~ x.fall + I(x.fall^2)))
    poly_mod <- lm(y.fall ~ x.fall + I(x.fall^2))     
    pred.y.fall <- predict(lm(y.fall ~ x.fall + I(x.fall^2)))
     fall.Jul.coeffs <- coefficients(Jul_summ_fall)
     fall.poly.coeffs <- coefficients(poly_summ_fall)
     

    pressstat_Jul <- PRESS(Jul_mod)
    pressstat_Jul$stat
    pressstat_noJul <- PRESS(noJul_mod)
    pressstat_noJul$stat
    pressstat_poly <- PRESS(poly_mod)
    pressstat_poly$stat
     
    metrics_out[1,7] <- pressstat_Jul$stat
    metrics_out[1,8] <- pressstat_Jul$P.square
    metrics_out[1,9] <- pressstat_noJul$stat
    metrics_out[1,10] <- pressstat_noJul$P.square
    metrics_out[1,11] <- pressstat_poly$stat
    metrics_out[1,12] <- pressstat_poly$P.square
     
     fall.resids <- residuals(Jul_summ_fall)
     write.table (x=fall.resids,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/resids_fall_poly_Mn.csv",sep = ",", col.names=F)       
     
     pred.out <- matrix(nrow=length(y.fall), ncol=3)
     pred.out[,1] <- y.fall
     pred.out[,2] <- pred.y.fall  
     pred.out[,3] <- as.character(NoNA.xyz$SiteName[maxrow:nrow(NoNA.xyz)])
     write.table (x=pred.out,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/Pred_fall_poly_Mn.csv",sep = ",", col.names=F)  
     remove(pred.out)
            
     
     coeffs_out[1,1] <- sp.Jul.coeffs[1,1]
     coeffs_out[1,2] <- sp.Jul.coeffs[2,1]
     coeffs_out[1,3] <- sp.Jul.coeffs[3,1]
     coeffs_out[1,4] <- fall.Jul.coeffs[1,1]
     coeffs_out[1,5] <- fall.Jul.coeffs[2,1]
     coeffs_out[1,6] <- fall.Jul.coeffs[3,1]
     coeffs_out[1,7] <- sp.poly.coeffs[1,1]
     coeffs_out[1,8] <- sp.poly.coeffs[2,1]
     coeffs_out[1,9] <- fall.poly.coeffs[1,1]
     coeffs_out[1,10] <- fall.poly.coeffs[2,1]
     coeffs_out[1,11] <- sp.poly.coeffs[3,1]
     coeffs_out[1,12] <- fall.poly.coeffs[3,1]
     
     stats_out[1,1] <- Jul_summ_sp$sigma
     stats_out[1,2] <- Jul_summ_sp$adj.r.squared
     stats_out[1,3] <- Jul_summ_fall$sigma
     stats_out[1,4] <- Jul_summ_fall$adj.r.squared
     stats_out[1,5] <- poly_summ_sp$sigma
     stats_out[1,6] <- poly_summ_sp$adj.r.squared
     stats_out[1,7] <- poly_summ_fall$sigma
     stats_out[1,8] <- poly_summ_fall$adj.r.squared

     
    

stats_out[,1:8] <- round(stats_out[,1:8], digits = 2)

write.table(x=stats_out, append=F,row.names=F, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/2012_mod_stats_Mn.csv", sep = ",", col.names=T)

write.table(x=coeffs_out, append=F,row.names=F, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/2012_mod_coeffs_Mn.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=F, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/2012_mod_metrics_Mn.csv", sep = ",", col.names=T)

write.table(x=LST.Log.out, append=F,row.names=F, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/2012_km_model_data_Mn.csv", sep = ",", col.names=T)

rm(list = ls())

########################################################################################################
# This part applies the model coefficients to the LST to generate predicted daily stream temps
# There are separate functions for using the LST-only coefficients and the LST + Julian day model coefficients
########################################################################################################

#first section has code for basin-wide model application

setwd("D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling")
LST.in <- read.dbf("D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/LST12_Ent_RCA.dbf")
colnames(LST.in)<-gsub("X", "", colnames(LST.in))
#newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    #specifies those dates to be numeric in order to drop the leading zeros
newnamesnum <- as.numeric(colnames(LST.in)[1:46])
colnames(LST.in)[1:46] <- newnamesnum;

coeffs.in <- read.csv("2012_mod_coeffs_Mn.csv")
LogPred.out <- LST.in
LogPred.out[,] <- 0


	for (i in 1:448)	#basin-wide model application - change for actual number of RCAs
            	{
            	x <- unlist(LST.in[i,])
            	maxrow <- as.numeric(which.max(x[1:46]))
            	midrow <- maxrow + 1
            
            	x[1:maxrow] <- x[1:maxrow] * coeffs.in$sp_poly_bLST[1] + x[1:maxrow]^2 * coeffs.in$sp_poly_bLST2[1] + coeffs.in$sp_poly_Int[1]
            	x[midrow:46] <- x[midrow:46] * coeffs.in$fall_Jul_bLST[1] + i * coeffs.in$fall_wJul_bJul[1]  + coeffs.in$fall_poly_Int[1] 
            	
            	LogPred.out[i,] <- x
		}

#for one with Julian day included, one with poly:

  for (i in 1:448)
    {
      x <- unlist(LST.in[i,])
      maxrow <- 29 #as.numeric(which.max(x[1:46]))
      midrow <- maxrow + 1
      day <- as.numeric(colnames(LST.in[midrow]))
    
      j <- 1
      for (l in 1:maxrow)
      {x[l] <- x[l] * coeffs.in$sp_poly_bLST[1] + x[l]^2 * coeffs.in$sp_poly_bLST2[1] + coeffs.in$sp_poly_Int[1]
       j <- j + 8}
      k <- day
      for (l in midrow:46)     
      {x[l] <- x[l] * coeffs.in$fall_Jul_bLST[1] + k * coeffs.in$fall_Jul_bJul[1] + coeffs.in$fall_Jul_Int[1]
       k <- k + 8}
      LogPred.out[i,] <- x  
    }

#for both wJul models:

  for (i in 1:448)
    {
        x <- unlist(LST.in[i,])
        maxrow <- 29 #as.numeric(which.max(x[1:46])) #either specify or let be dynamic
        midrow <- maxrow + 1
        day <- as.numeric(colnames(LST.in[midrow]))
        
        j <- 1
        for (l in 1:maxrow)
        {x[l] <- x[l] * coeffs.in$sp_Jul_bLST[1] + j * coeffs.in$sp_Jul_bJul[1] + coeffs.in$sp_Jul_Int[1]
         j <- j + 8}
        k <- day
        for (l in midrow:46)     
        {x[l] <- x[l] * coeffs.in$fall_Jul_bLST[1] + k * coeffs.in$fall_Jul_bJul[1] + coeffs.in$fall_Jul_Int[1]
         k <- k + 8}
        LogPred.out[i,] <- x  
    }

LogPred.out <- as.data.frame(LogPred.out)
LogPred.out[LogPred.out<0] = 0.0
#LogPred.out[LogPred.out[,1:46]>50] = "NA"
#colnames(LogPred.out)[2:47] <- day.names[1:46]
write.dbf(LogPred.out, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/23Aug14modeling/predt2012_Ent_8D_Mn.dbf")

#this section has a function for more complicated basin partitioning model application

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

write.dbf(merged.LogPred, file = "D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/2Feb14modelling/predt2013_JD_8D_Mn.dbf")



rm(list = ls())


#################################################################
#Summary stats - reads in a years worth of 8-day predictions and summarizes for various annual and seEntnal periods
##################################################################

meanTemp <- read.dbf("C:/Dropbox/work/research/CHaMP/CHaMP_data/LST/predt2012_JD_Mn.dbf")
maxTemp <- read.dbf("C:/Dropbox/work/research/CHaMP/CHaMP_data/LST/predt2012_JD_Mx.dbf")
minTemp <- read.dbf("C:/Dropbox/work/research/CHaMP/CHaMP_data/LST/predt2012_JD_Min.dbf")

pred_summary <- data.frame(meanTemp[,1])


for (i in 1:5532)
     {x <- unlist(maxTemp[i,2:35])
     xstat <- max(x, na.rm = TRUE)
     pred_summary$AnnMxMx[i] <- xstat
     }

for (i in 1:5532)
     {x <- unlist(maxTemp[i,25:35])
     xstat <- mean(x, na.rm = TRUE)
     pred_summary$JulSeptMnMx[i] <- xstat
     }

for (i in 1:5532)
     {x <- unlist(meanTemp[i,25:35])
     xstat <- mean(x, na.rm = TRUE)
     pred_summary$JulSeptMnMn[i] <- xstat
     }

for (i in 1:5532)
     {x <- unlist(meanTemp[i,2:35])
     xstat <- mean(x, na.rm = TRUE)
     pred_summary$AnnMnMn[i] <- xstat
     }

for (i in 1:5532)
     {x <- unlist(minTemp[i,c(2:21,30:35)])
     xstat <- mean(x, na.rm = TRUE)
     pred_summary$OctMayMnMin[i] <- xstat
     }

for (i in 1:5532)
     {x <- unlist(meanTemp[i,9:25])
     xstat <- mean(x, na.rm = TRUE)
     pred_summary$MarJunMnMn[i] <- xstat
     }
     
pred_summary$ObjID_JD_s <- maxTemp$ObjID_JD_s
pred_summary$OBJECTID <- maxTemp$OBJECTID

write.table(x=pred_summary, append=F,row.names=F, file = "C:/Dropbox/work/research/CHaMP/CHaMP_data/LST/JD_2012_summary.csv", sep = ",", col.names=T)

#####################################

Fpar<- read.dbf("C:/Dropbox/work/research/CHaMP/CHaMP_data/2012/JohnDay/Solar/kmpt_JD_Fpar_2012.dbf")
Fpar[Fpar>249] = NA
LST.in[LST.in<0.1] = NA






####################################################################################
#This part reads in the training dataset for a year and generates RMSEs for site-specific models
#####################################################################################

RcaID.in <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/RCA_DCE_unique.csv")
Log.in <- read.csv("JohnDay_8Day_2013.csv")
rca_elev <- read.csv("C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/RCA_zone_min_elev.csv")

pathname <-"C:/Dropbox/work/research/Steelhead/data_analysis/LST/predictions/site_regressions/2009/1Oct12/"
RcaID.in <- merge(RcaID.in, rca_elev, by.x = "RASTERVALU", by.y = "RASTERVALU", all.x = TRUE, all.y = FALSE)

setwd("D:/Dropbox/work/research/CHaMP/CHaMP_data/JohnDay/31Jan13modelling/")
LST.in <- read.dbf("D:/Dropbox/work/research/CHaMP/CHaMP_data/LST/LST12_JD_RCA_smoothed.dbf")

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
     
    write.table (x=pred.out,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/18Aug14modeling/Y_predY_2013_noJul_train.csv",sep = ",", col.names=F)
    write.table (x=mergeLogLST,append=T,row.names=F,file="D:/Dropbox/work/research/CHaMP/CHaMP_data/Entiat/18Aug14modeling/2013_noJul_train_preds.csv",sep = ",", col.names=F)
    
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


