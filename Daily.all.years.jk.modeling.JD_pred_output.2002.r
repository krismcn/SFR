############################################################################################################
# This set of R scripts runs a by-site LOO jk on the daily JD data for model validation
# Created: 14 July 2015

          

##############################################################################################


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

 
  
########
  setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/JohnDay/")
  
  ID.in <- read.csv("OBJID_RCAID_RCA_shape.csv")


  
  rca_elev <- read.csv("D:/OneDrive/work/research/Steelhead/spatial_data/rca_elev.csv")
  
   LST.in <- read.dbf("D:/OneDrive/work/research/Steelhead/spatial_data/Linear_interp/RCAzone_LST_02.dbf")
  colnames(LST.in)<-gsub("X", "", colnames(LST.in))
  #newnamesnum <- as.numeric(substring(colnames(LST.in)[1:46],3,5));    #specifies those dates to be numeric in order to drop the leading zeros
  newnamesnum <- as.numeric(colnames(LST.in)[1:365]);
  colnames(LST.in)[1:365] <- newnamesnum;
  #LST.in[,1:365] <- LST.in[,1:365]-273.15
  
   HUCs.in <- read.csv("D:/Dropbox/work/research/Steelhead/spatial_data/coverages/RCA_HUCs.csv")
  
  
  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/")
  
  
  ID <- merge(ID.in, HUCs.in, by.x = "OBJECTID", by.y = "OBJECTID", all.x = TRUE, all.y = FALSE)
  ID.in <- ID
  remove(ID)

  LST <- merge(LST.in, ID.in, by.x = "OBJECTID", by.y = "OBJECTID", )

  q <- "2002"

  data.sp <- read.csv("All_model_data_sp_2002.csv", header=F, stringsAsFactors = F)
  colnames(data.sp) <- c("y", "x", "z", "e", "HUC", "SiteName")
  

  data.fall <- read.csv("All_model_data_fall_2002.csv", header=F, stringsAsFactors = F)
  colnames(data.fall) <- c("y", "x", "z", "e", "HUC", "SiteName")
  

  basin <- "JD"
  longBasin <- "JohnDay"
  yearPath <- "2002"
  
  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/"
  
  setwd(paste0(mainPath, longBasin, "/modelPreds/", yearPath))

  
###########################################
# HUC models
###########################################

  coeffs_out <- data.frame(Int=numeric(0), bLST=numeric(0), bLST2=numeric(0), bJul=numeric(0), bElev=numeric(0), HUC=numeric(0))
  metrics_out <- data.frame(PRESS=numeric(0), p2=numeric(0), r2=numeric(0), RMSEP=numeric(0), RMSE=numeric(0), HUC=numeric(0))
  
  write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mn.csv"), sep = ",", col.names=NA)
  write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_Mn.csv"), sep = ",", col.names=NA)
  

    for (j in 1:4)
    {
      
      HUC.sp <- data.sp[data.sp$HUC == j,]
      
        y <- HUC.sp$y
        x <- HUC.sp$x
        z <- HUC.sp$z
        e <- HUC.sp$e
        plot(x, y, main = j)
        plot(z, y, main = j)
        mod <- lm(y ~ x + I(x^2) + z + e)
        sum_mod <- summary(mod)
        coeffs <- as.matrix(coefficients(mod))
        pred.y <- predict(mod)
        pred.y[pred.y<0] = 0.0
        plot(y, pred.y, main = j)
  
  
        pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
        RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
        RSS <- sum(pressstat_sum$residuals^2)
        TSS <- sum((y - mean(y))^2)
        MSS <- sum((pred.y - mean(y))^2)
        p2 <- MSS/TSS
        
        pred.out <- matrix(nrow=length(pred.y), ncol=6)
        pred.out[,1] <- y
        pred.out[,2] <- pred.y
        pred.out[,3] <- z
        pred.out[,4] <- "spring"
        pred.out[,5] <- yearPath
        pred.out[,6] <- j
        colnames(pred.out) <- c("Y", "PredY", "JulDay", "Season", "Year", "HUC")
        write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall_all.csv"), sep = ",", col.names=F)  
  
  
  
        coeffs_out[1,1] <- coeffs[1,1]
        coeffs_out[1,2] <- coeffs[2,1]
        coeffs_out[1,3] <- coeffs[3,1]
        coeffs_out[1,4] <- coeffs[4,1]
        coeffs_out[1,5] <- coeffs[5,1]
        coeffs_out[1,6] <- j
        
        metrics_out[1,1] <- pressstat_sum$stat
        metrics_out[1,2] <- p2
        metrics_out[1,3] <- sum_mod$adj.r.squared
        metrics_out[1,4] <- RMSEP
        metrics_out[1,5] <- sum_mod$sigma
        metrics_out[1,6] <- j
  
        HUC.fall <- data.fall[data.fall$HUC == j,]
        y <- HUC.fall$y
        x <- HUC.fall$x
        z <- HUC.fall$z
        e <- HUC.fall$e
        plot(x, y, main = j)
        plot(z, y, main = j)
        mod <- lm(y ~ x + I(x^2) + z + e)
        sum_mod <- summary(mod)
        coeffs <- as.matrix(coefficients(mod))
        pred.y <- predict(mod)
        pred.y[pred.y<0] = 0.0
        plot(y, pred.y, main = j)
  
  
        pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
        RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))
        RSS <- sum(pressstat_sum$residuals^2)
        TSS <- sum((y - mean(y))^2)
        MSS <- sum((pred.y - mean(y))^2)
        p2 <- MSS/TSS
  
        pred.out <- matrix(nrow=length(pred.y), ncol=6)
        pred.out[,1] <- y
        pred.out[,2] <- pred.y
        pred.out[,3] <- z
        pred.out[,4] <- "fall"
        pred.out[,5] <- yearPath
        pred.out[,6] <- j
      
        write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall_all.csv"),sep = ",", col.names=F)  
  
        
  
        coeffs_out[2,1] <- coeffs[1,1]
        coeffs_out[2,2] <- coeffs[2,1]
        coeffs_out[2,3] <- coeffs[3,1]
        coeffs_out[2,4] <- coeffs[4,1]
        coeffs_out[2,5] <- coeffs[5,1]
        coeffs_out[2,6] <- j
        
        metrics_out[2,1] <- pressstat_sum$stat
        metrics_out[2,2] <- p2
        metrics_out[2,3] <- sum_mod$adj.r.squared
        metrics_out[2,4] <- RMSEP
        metrics_out[2,5] <- sum_mod$sigma
        metrics_out[2,6] <- j
        
        rownames(metrics_out) <- c("Spring", "Fall")
        rownames(coeffs_out) <- c("Spring", "Fall")
        write.table(x=coeffs_out, append=T,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mn.csv"), sep = ",", col.names=F)
        
        write.table(x=metrics_out, append=T,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_Mn.csv"), sep = ",", col.names=F)
  
    }


  # model predictions
  ######################################################
  library(timeSeries)
  library(lattice)
  library(foreign)
  library(doBy)
  library(qpcR)
  library(pls)
  library(boot)
  library(Hmisc)
  library(zoo)
  
  #rca_elev <- read.csv("D:/OneDrive/work/research/Steelhead/spatial_data/rca_elev.csv")
  
  
  LST.elev <- merge(LST, rca_elev, by.x = "RCA_ID", by.y = "RCA_ID")
  
  setwd(paste0(mainPath, longBasin, "/modelPreds/", yearPath))

  coeffs.in <- read.csv(paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_Mn.csv"), stringsAsFactors = FALSE)
 
  
#   for (i in 1:length(rcas))  
#   {
#     x <- unlist(LST.elev[i,])
#     maxrow <- as.numeric(which.max(x[3:367]))
#     midrow <- maxrow + 1
#     j <- 1
#     day <- as.numeric(colnames(LST.elev)[maxrow])
#     
#     
#     for (l in 1:maxrow)
#     {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[371] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
#      j <- j + 8}
#     k <- midrow
#     for (l in midrow:367)     
#     {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[370] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
#      k <- k + 8}
#     LogPred.out[i,1:12] <- x [1:12] 
#   }

  
  
  #this section has a function for more complicated basin partitioning model application
  
  #HUCs.in <- read.csv("D:/Dropbox/work/research/Steelhead/spatial_data/coverages/RCA_HUCs.csv")
  #LST.HUCs <- merge(LST.in, HUCs.in, by.x = "OBJECTID", by.y = "OBJECTID", all.x = F, all.y = T)
#   LST.HUCs.out <- LST.elev
#   LST.HUCs.out[,] <- 0
#   
  
  myFunc <- function(x) #HUC based model application
  {
    maxrow <- as.numeric(which.max(x[3:367]))
    midrow <- maxrow + 1
    j <- 1
    
    if (x[368] == 1)
    {
      for (l in 3:maxrow)
      {x[l] <- x[l] * coeffs.in$bLST[1] + j * coeffs.in$bJul[1] +  x[370] * coeffs.in$bElev[1] + x[l]^2 * coeffs.in$bLST2[1] + coeffs.in$Int[1]
       j <- j + 1}
      k <- midrow
      for (l in midrow:367)     
      {x[l] <- x[l] * coeffs.in$bLST[2] + k * coeffs.in$bJul[2] +  x[370] * coeffs.in$bElev[2] + x[l]^2 * coeffs.in$bLST2[2] + coeffs.in$Int[2]
       k <- k + 1}
    } else
      if (x[368] == 2)
      {
        for (l in 3:maxrow)
        {x[l] <- x[l] * coeffs.in$bLST[3] + j * coeffs.in$bJul[3] +  x[370] * coeffs.in$bElev[3] + x[l]^2 * coeffs.in$bLST2[3] + coeffs.in$Int[3]
         j <- j + 1}
        k <- midrow
        for (l in midrow:367)     
        {x[l] <- x[l] * coeffs.in$bLST[4] + k * coeffs.in$bJul[4] +  x[370] * coeffs.in$bElev[4] + x[l]^2 * coeffs.in$bLST2[4] + coeffs.in$Int[4]
         k <- k + 1}
      } else
        if (x[368] == 3)
        {
          for (l in 3:maxrow)
          {x[l] <- x[l] * coeffs.in$bLST[5] + j * coeffs.in$bJul[5] +  x[370] * coeffs.in$bElev[5] + x[l]^2 * coeffs.in$bLST2[5] + coeffs.in$Int[5]
           j <- j + 1}
          k <- midrow
          for (l in midrow:367)     
          {x[l] <- x[l] * coeffs.in$bLST[6] + k * coeffs.in$bJul[6] +  x[370] * coeffs.in$bElev[6] + x[l]^2 * coeffs.in$bLST2[6] + coeffs.in$Int[6]
           k <- k + 1}
        } else
          if (x[368] == 4)
          {
            for (l in 3:maxrow)
            {x[l] <- x[l] * coeffs.in$bLST[7] + j * coeffs.in$bJul[7] +  x[370] * coeffs.in$bElev[7] + x[l]^2 * coeffs.in$bLST2[7] + coeffs.in$Int[7]
             j <- j + 1}
            k <- midrow
            for (l in midrow:367)     
            {x[l] <- x[l] * coeffs.in$bLST[8] + k * coeffs.in$bJul[8] +  x[370] * coeffs.in$bElev[8] + x[l]^2 * coeffs.in$bLST2[8] + coeffs.in$Int[8]
             k <- k + 1}
            #fill <- rollmean(x[3:367], 5, align = "left")
            #x[180:210] <- fill[180:210]
          }
    x
  }

  LogPred.out <- apply(LST.elev, 1, myFunc)
  tLogPred.out <- as.data.frame(t(LogPred.out))
  tLogPred.out[tLogPred.out<0] = 0.0
  
  tLogPred.out[,3:367] <- round(tLogPred.out[,3:367], digits = 2)
  
  plot(1:365, tLogPred.out[3000,3:367])#HUC 4
  points(1:365, tLogPred.out[100,3:367], pch=16, col="blue")#HUC 2
  points(1:365, tLogPred.out[2000,3:367], pch=16, col="lightblue")#HUC 3
  points(1:365, tLogPred.out[5000,3:367], pch=16, col="purple")#HUC 1


  write.dbf(tLogPred.out, file = "predt2002_JD_daily_Mn.dbf")
  
#####################################
# spring/fall basin-wide
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
  plot(x, y, main = "all spring")
  plot(z, y, main = "all spring")
  
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
  write.table (x=pred.out,append=F,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall_basin.csv"), sep = ",", col.names=F)  
  
  
  plot(y, pred.y, main = "all spring")
  
  
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
  plot(x, y, main = "All fall")
  plot(z, y, main = "All fall")
  mod <- lm(y ~ x + I(x^2) + z + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))
  pred.y <- predict(mod)
  pred.y[pred.y<0] = 0.0
  plot(y, pred.y, main = "All fall")
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
  
  write.table (x=pred.out,append=T,row.names=F,file=paste0("jk_pred_v_y_", basin, "_", yearPath, "_sp_fall_all.csv"),sep = ",", col.names=F)  
  
  
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
  
  write.table(x=coeffs_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_coeffs_basin_Mn.csv"), sep = ",", col.names=T)
  
  write.table(x=metrics_out, append=F,row.names=T, file = paste0("All_data_", basin, "_", yearPath, "_mod_metrics_basin_Mn.csv"), sep = ",", col.names=T)
