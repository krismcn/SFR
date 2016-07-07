############################################################################################################
# This set of R scripts to process NorWeST data and runs a by-site LOO jk on the daily summer only JD data for model validation and summer model comparison
# Created: 24 July 2015

          

##############################################################################################

library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/")

sites <- read.csv("Clearwater_CBS_sites.csv", header=T)
data <- read.csv("Clearwater_CBS_logger_data.csv", header=T)
sites$Orig_ID <- as.character(sites$Orig_ID)
data$OrigID <- as.character(data$OrigID)
all_data <- merge(data, sites, by.x= "OrigID", by.y = "Orig_ID", all.x = TRUE, all.y = FALSE)
data_2011 <- all_data[all_data$Year == 2011,]
data_2012 <- all_data[all_data$Year == 2012,]
data_2013 <- all_data[all_data$Year == 2013,]

all_data$OrigID <- as.character(all_data$OrigID)
ID.2011 <- as.data.frame(table(data_2011$OrigID)) 
colnames(NoNA.xyz) <- c("y", "x", "z", "e", "HUC", "SiteName")

  summer <- subset(NoNA.xyz, z > 181 & z < 258)
  SiteID <- as.matrix(unique(summer$SiteName))
  SiteID <- as.matrix(as.character(SiteID))
  ID.table <- as.data.frame(table(summer$SiteName))  
  SiteID.sum <- as.matrix(ID.table[ID.table$Freq > 40, "Var1"])
 
setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/summer_all_years")  
   
  pred.out <- matrix(nrow=0, ncol = 6)
  colnames(pred.out) <- c("Y", "PredY", "JulDay", "SiteName", "season", "Year")
  write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_JD_2004.csv",sep = ",", col.names=T)
  stats.out <- matrix(nrow=0, ncol = 10)
  colnames(stats.out) <- c("SiteName", "Year", "RMSE", "r2", "Int", "x", "x2", "z", "e", "RMSEP")
  write.table (x=stats.out, append=F, row.names=F, file="jk_pred_metrics_JD_2004.csv", sep=",", col.names=T)

  j <- 0
  
    while(j <= 100)
      {        
        sites <- sample(SiteID.sum, 16, replace = FALSE, prob=NULL) #change sample size with each data subset
        
        data <- summer[summer$SiteName %in% sites,]
        test.data <- summer[summer$SiteName %nin% sites,]
        
        y <- data$y
        x <- data$x
        z <- data$z
        e <- data$e
        mod <- lm(y ~ x + I(x^2) + z + e)
        sum_mod <- summary(mod)
        coeffs <- as.matrix(coefficients(mod))
        pred.y <- predict(mod, newdata = test.data)
        pred.y[pred.y< -0.5] = -0.5
        RMSEP <- sqrt(mean((test.data$y - pred.y)^2))
        
        
        pred.out <- matrix(nrow=length(pred.y), ncol=6)
        pred.out[,1] <- test.data$y
        pred.out[,2] <- pred.y
        pred.out[,3] <- test.data$z
        pred.out[,4] <- as.character(test.data$SiteName)
        pred.out[,5] <- "summer"
        pred.out[,6] <- "2004"
      
        write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD_2004.csv",sep = ",", col.names=F)  
        
        stats.out <- matrix(nrow=1, ncol=10)
        stats.out[1,1] <- as.character(test.data$SiteName)[1]
        stats.out[1,2] <- "2004"
        stats.out[1,3] <- sum_mod$sigma
        stats.out[1,4] <- sum_mod$adj.r.squared
        stats.out[1,5] <- coeffs[1,1]
        stats.out[1,6] <- coeffs[2,1]
        stats.out[1,7] <- coeffs[3,1]
        stats.out[1,8] <- coeffs[4,1]
        stats.out[1,9] <- coeffs[5,1]
        stats.out[1,10] <- RMSEP
        
        write.table (x=stats.out, append=T, row.names=F, file="jk_pred_metrics_JD_2004.csv", sep=",", col.names=F)
        
        j <- j+1
      }


######################################################
# no Julday
#####

pred.out <- matrix(nrow=0, ncol = 6)
colnames(pred.out) <- c("Y", "PredY", "JulDay", "SiteName", "season", "Year")
write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_JD_2004_noJul.csv",sep = ",", col.names=T)
stats.out <- matrix(nrow=0, ncol = 9)
colnames(stats.out) <- c("SiteName", "Year", "RMSE", "r2", "Int", "x", "x2", "e", "RMSEP")
write.table (x=stats.out, append=F, row.names=F, file="jk_pred_metrics_JD_2004_noJul.csv", sep=",", col.names=T)

j <- 0

  while(j <= 100)
    {        
      sites <- sample(SiteID.sum, 16, replace = FALSE, prob=NULL) #change sample size with each data subset
      
      data <- summer[summer$SiteName %in% sites,]
      test.data <- summer[summer$SiteName %nin% sites,]
      
      y <- data$y
      x <- data$x
      z <- data$z
      e <- data$e
      mod <- lm(y ~ x + I(x^2) + e)
      sum_mod <- summary(mod)
      coeffs <- as.matrix(coefficients(mod))
      pred.y <- predict(mod, newdata = test.data)
      RMSEP <- sqrt(mean((test.data$y - pred.y)^2))
      pred.y[pred.y<0] = 0.0
      
      pred.out <- matrix(nrow=length(pred.y), ncol=6)
      pred.out[,1] <- test.data$y
      pred.out[,2] <- pred.y
      pred.out[,3] <- test.data$z
      pred.out[,4] <- as.character(test.data$SiteName)
      pred.out[,5] <- "summer"
      pred.out[,6] <- "2004"
      
      write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD_2004_noJul.csv",sep = ",", col.names=F)  
      
      stats.out <- matrix(nrow=1, ncol=10)
      stats.out[1,1] <- as.character(test.data$SiteName)[1]
      stats.out[1,2] <- "2004"
      stats.out[1,3] <- sum_mod$sigma
      stats.out[1,4] <- sum_mod$adj.r.squared
      stats.out[1,5] <- coeffs[1,1]
      stats.out[1,6] <- coeffs[2,1]
      stats.out[1,7] <- coeffs[3,1]
      stats.out[1,8] <- coeffs[4,1]
      stats.out[1,9] <- RMSEP
      
      write.table (x=stats.out, append=T, row.names=F, file="jk_pred_metrics_JD_2004_noJul.csv", sep=",", col.names=F)
      
      j <- j+1
    }

###########################################################################################
# Full year models




pred.out <- matrix(nrow=0, ncol = 6)
colnames(pred.out) <- c("Y", "PredY", "JulDay", "SiteName", "season", "Year")
write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_JD_2004_yr.csv",sep = ",", col.names=T)
stats.out <- matrix(nrow=0, ncol = 9)
colnames(stats.out) <- c("SiteName", "Year", "RMSE", "r2", "Int", "x", "x2", "e", "RMSEP")
write.table (x=stats.out, append=F, row.names=F, file="jk_pred_metrics_JD_2004_yr.csv", sep=",", col.names=T)

  NoNA.xyz <- orderBy(~z, NoNA.xyz)
  maxrow <- which.max(NoNA.xyz$y)
  NoNA.xyz[maxrow,"z"]
  
  spring.data <- NoNA.xyz[1:maxrow,c(1,2,3,4,6)]
  fall.data <- NoNA.xyz[maxrow:nrow(NoNA.xyz),c(1,2,3,4,6)]
  spring.data <- orderBy(~SiteName, spring.data)
  fall.data <- orderBy(~SiteName, fall.data) 
  
  spring.data$SiteName <- as.character(spring.data$SiteName)
  sp.table <- as.data.frame(table(spring.data$SiteName))  
  SiteID.sp <- as.matrix(sp.table[sp.table$Freq > 40, "Var1"])
  
  fall.data$SiteName <- as.character(fall.data$SiteName)
  fall.table <- as.data.frame(table(fall.data$SiteName))  
  SiteID.fall <- as.matrix(fall.table[fall.table$Freq > 40, "Var1"])
  
  
 
  j <- 0

  while(j <= 100)
    {        
     
      sites.sp <- sample(SiteID.sp, 9, replace = FALSE, prob=NULL) #change sample size with each data subset
      
      data.sp <- spring.data[spring.data$SiteName %in% sites.sp,]
      test.data.sp <- spring.data[spring.data$SiteName %nin% sites.sp,]
      
      y <- data.sp$y
      x <- data.sp$x
      z <- data.sp$z
      e <- data.sp$e
      mod <- lm(y ~ x + I(x^2) + z + e)
      sum_mod <- summary(mod)
      coeffs <- as.matrix(coefficients(mod))
      pred.y <- predict(mod, newdata = test.data.sp)
      RMSEP <- sqrt(mean((test.data.sp$y - pred.y)^2))
      pred.y[pred.y<0] = 0.0
      
      pred.out <- matrix(nrow=length(pred.y), ncol=6)
      pred.out[,1] <- test.data.sp$y
      pred.out[,2] <- pred.y
      pred.out[,3] <- test.data.sp$z
      pred.out[,4] <- as.character(test.data.sp$SiteName)
      pred.out[,5] <- "spring"
      pred.out[,6] <- "2004"

      write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD_2004_yr.csv",sep = ",", col.names=F)  
      
      stats.out <- matrix(nrow=1, ncol=10)
      stats.out[1,1] <- as.character(test.data.sp$SiteName)[1]
      stats.out[1,2] <- "2004"
      stats.out[1,3] <- sum_mod$sigma
      stats.out[1,4] <- sum_mod$adj.r.squared
      stats.out[1,5] <- coeffs[1,1]
      stats.out[1,6] <- coeffs[2,1]
      stats.out[1,7] <- coeffs[3,1]
      stats.out[1,8] <- coeffs[4,1]
      stats.out[1,9] <- RMSEP
      
      write.table (x=stats.out, append=T, row.names=F, file="jk_pred_metrics_JD_2004_yr.csv", sep=",", col.names=F)
      
      
      sites.fall <- sample(SiteID.fall, 16, replace = FALSE, prob=NULL) #change sample size with each data subset
      
      data.fall <- spring.data[spring.data$SiteName %in% sites.fall,]
      test.data.fall <- spring.data[spring.data$SiteName %nin% sites.fall,]
      
      y <- data.fall$y
      x <- data.fall$x
      z <- data.fall$z
      e <- data.fall$e
      mod <- lm(y ~ x + I(x^2) + z + e)
      sum_mod <- summary(mod)
      coeffs <- as.matrix(coefficients(mod))
      pred.y <- predict(mod, newdata = test.data.fall)
      RMSEP <- sqrt(mean((test.data.fall$y - pred.y)^2))
      pred.y[pred.y<0] = 0.0
      
      pred.out <- matrix(nrow=length(pred.y), ncol=6)
      pred.out[,1] <- test.data.fall$y
      pred.out[,2] <- pred.y
      pred.out[,3] <- test.data.fall$z
      pred.out[,4] <- as.character(test.data.fall$SiteName)
      pred.out[,5] <- "fall"
      pred.out[,6] <- "2004"
      
      write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD_2004_yr.csv",sep = ",", col.names=F)  
      
      stats.out <- matrix(nrow=1, ncol=10)
      stats.out[1,1] <- as.character(test.data.fall$SiteName)[1]
      stats.out[1,2] <- "2004"
      stats.out[1,3] <- sum_mod$sigma
      stats.out[1,4] <- sum_mod$adj.r.squared
      stats.out[1,5] <- coeffs[1,1]
      stats.out[1,6] <- coeffs[2,1]
      stats.out[1,7] <- coeffs[3,1]
      stats.out[1,8] <- coeffs[4,1]
      stats.out[1,9] <- RMSEP
      
      write.table (x=stats.out, append=T, row.names=F, file="jk_pred_metrics_JD_2004_yr.csv", sep=",", col.names=F)
      
      j <- j+1
    }
  

