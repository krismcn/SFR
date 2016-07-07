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

 
  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/jk_all_years")  
  
  pred.out <- matrix(nrow=0, ncol = 6)
  colnames(pred.out) <- c("Y", "PredY", "SiteName", "season", "Year", "SE")
  write.table (x=pred.out,append=F,row.names=F,file="jk_pred_v_y_JD.csv",sep = ",", col.names=T)
  stats.out <- matrix(nrow=0, ncol = 10)
  colnames(stats.out) <- c("Year", "Season", "RMSE", "r2", "Int", "x", "x2", "z", "e", "RMSEP")
  write.table (x=stats.out, append=F, row.names=F, file="jk_pred_RMSEP_JD.csv", sep=",", col.names=T)

########

  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/")
  q <- "2009"

  spring.data <- read.csv("All_model_data_sp_2009.csv", header=F)
  colnames(spring.data) <- c("y", "x", "z", "e", "HUC", "SiteName")
  spring.data$SiteName <- as.character(spring.data$SiteName)

  sp.table <- as.data.frame(table(spring.data$SiteName))  
  SiteID.sp <- as.matrix(sp.table[sp.table$Freq > 10, "Var1"])
  site.num <- nrow(SiteID.sp)
 
  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/jk_all_years")

  j <- 0
  
    while(j <= 100)
      {        
        sites <- sample(SiteID.sp, site.num-1, replace = FALSE, prob=NULL) #change sample size with each data subset
        
        data <- spring.data[spring.data$SiteName %in% sites,]
        test.data <- spring.data[spring.data$SiteName %nin% sites,]
        
        y <- data$y
        x <- data$x
        z <- data$z
        e <- data$e
        mod <- lm(y ~ x + I(x^2) + z + e)
        sum_mod <- summary(mod)
        coeffs <- as.matrix(coefficients(mod))
        
        pred.y.spring <- predict(mod, newdata = test.data)
        pred.y.spring[pred.y.spring< -0.5] = -0.5
        RMSEP <- sqrt(mean((test.data$y - pred.y.spring)^2))
        SE <- (test.data$y - pred.y.spring)^2
        
        pred.out <- matrix(nrow=length(pred.y.spring), ncol=6)
        pred.out[,1] <- test.data$y
        pred.out[,2] <- pred.y.spring
        pred.out[,3] <- as.character(test.data$SiteName)
        pred.out[,4] <- "spring"
        pred.out[,5] <- q
        pred.out[,6] <- SE
      
        write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD.csv",sep = ",", col.names=F)  
        
        stats.out <- matrix(nrow=1, ncol=10) 
        stats.out[1,1] <- q
        stats.out[1,2] <- "spring"
        stats.out[1,3] <- sum_mod$sigma
        stats.out[1,4] <- sum_mod$adj.r.squared
        stats.out[1,5] <- coeffs[1,1]
        stats.out[1,6] <- coeffs[2,1]
        stats.out[1,7] <- coeffs[3,1]
        stats.out[1,8] <- coeffs[4,1]
        stats.out[1,9] <- coeffs[5,1]
        stats.out[1,10] <- RMSEP
        
        write.table (x=stats.out, append=T, row.names=F, file="jk_pred_RMSEP_JD.csv", sep=",", col.names=F)
        
        j <- j+1
      }

#############Fall

  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/")

  fall.data <- read.csv("All_model_data_fall_2009.csv", header=F)
  colnames(fall.data) <- c("y", "x", "z", "e", "HUC", "SiteName")
  fall.data$SiteName <- as.character(fall.data$SiteName)

  fall.table <- as.data.frame(table(fall.data$SiteName))  
  SiteID.fall <- as.matrix(fall.table[fall.table$Freq > 10, "Var1"])
  site.num <- nrow(SiteID.fall)

  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/jk_all_years")

    j <- 0
        
          while(j <= 100)
            {        
              sites <- sample(SiteID.fall, site.num-1, replace = FALSE, prob=NULL) #change sample size with each data subset
                
              data <- fall.data[fall.data$SiteName %in% sites,]
              test.data <- fall.data[fall.data$SiteName %nin% sites,]
                
              y <- data$y
              x <- data$x
              z <- data$z
              e <- data$e
              mod <- lm(y ~ x + I(x^2) + z + e)
              sum_mod <- summary(mod)
              coeffs <- as.matrix(coefficients(mod))
                
              pred.y.fall <- predict(mod, newdata = test.data)
              pred.y.fall[pred.y.fall< -0.5] = -0.5
              RMSEP <- sqrt(mean((test.data$y - pred.y.fall)^2))
              SE <- (test.data$y - pred.y.fall)^2
              
              pred.out <- matrix(nrow=length(pred.y.fall), ncol=6)
              pred.out[,1] <- test.data$y
              pred.out[,2] <- pred.y.fall
              pred.out[,3] <- as.character(test.data$SiteName)
              pred.out[,4] <- "fall"
              pred.out[,5] <- q
              pred.out[,6] <- SE
              
              write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD.csv",sep = ",", col.names=F)  
              
              stats.out <- matrix(nrow=1, ncol=10) 
              stats.out[1,1] <- q
              stats.out[1,2] <- "fall"
              stats.out[1,3] <- sum_mod$sigma
              stats.out[1,4] <- sum_mod$adj.r.squared
              stats.out[1,5] <- coeffs[1,1]
              stats.out[1,6] <- coeffs[2,1]
              stats.out[1,7] <- coeffs[3,1]
              stats.out[1,8] <- coeffs[4,1]
              stats.out[1,9] <- coeffs[5,1]
              stats.out[1,10] <- RMSEP
              
              write.table (x=stats.out, append=T, row.names=F, file="jk_pred_RMSEP_JD.csv", sep=",", col.names=F)
              
              j <- j+1
            }
          

#############All data all years

setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/")

all.data <- matrix(nrow=0, ncol = 8)
colnames(all.data) <- c("y", "x", "z", "e", "HUC", "SiteName", "Season", "Year")
write.table (x=all.data, append=F, row.names=F, file="all_train_data_all_years.csv", sep=",", col.names=T)

  for (i in 2000:2009)
    {
      spring.data <- read.csv(paste0("All_model_data_sp_", i, ".csv"), header=F)
      spring.data$Season <- "spring"
      spring.data$Year <- i
      fall.data <- read.csv(paste0("All_model_data_fall_", i, ".csv"), header=F)
      fall.data$Season <- "fall"
      fall.data$Year <- i
      write.table (x=spring.data, append=T, row.names=F, file="all_train_data_all_years.csv", sep=",", col.names=F)
      write.table (x=fall.data, append=T, row.names=F, file="all_train_data_all_years.csv", sep=",", col.names=F)
    }

all.data <- read.csv("All_train_data_all_years.csv", header=T)
spring.data <- all.data[all.data$Season == "spring",]
fall.data <- all.data[all.data$Season == "fall",]

    #data <- fall.data[fall.data$SiteName %in% sites,]
    #test.data <- fall.data[fall.data$SiteName %nin% sites,]
    
    y <- all.data$y
    x <- all.data$x
    z <- all.data$z
    e <- all.data$e
    year <- all.data$Year

    mod <- lm(y ~ x + I(x^2) + z + e)
    sum_mod <- summary(mod)
    coeffs <- as.matrix(coefficients(mod))
    
    all_mod <- lm(y ~ x + I(x^2) + z + e + x*z + x*e + e*z)
    noJul_mod <- lm(y ~ x)
    noJul_summ <- summary(noJul_mod)
    
    
    Jul_mod <- lm(y ~ x + z)
    Jul_summ <- summary(Jul_mod)
    poly_mod <- lm(y ~ x + I(x^2))
    poly_summ <- summary(poly_mod)
    elev_mod <- lm(y ~ x + z + e)
    elev_summ <- summary(elev_mod)     
    
    pred.y.spring.poly <- predict(lm(y.spring ~ x.spring + I(x.spring^2)))
    sp.poly.coeffs <- coefficients(poly_summ_sp)
    sp.Jul.coeffs <- coefficients(Jul_summ_sp)
    sp.resids <- residuals(poly_summ_sp)
    sp.elev.coeffs <- coefficients(elev_summ_sp)
    
    pressstat_all <- PRESS(mod)
    pressstat_Jul <- PRESS(Jul_mod)
    pressstat_noJul <- PRESS(noJul_mod)
    pressstat_poly <- PRESS(poly_mod)
    pressstat_elev <- PRESS(elev_mod_sp)
    
    pressstat_all$stat
    pressstat_Jul$stat
    pressstat_noJul$stat
    pressstat_poly$stat
    pressstat_elev$stat

    AIC(mod)
    AIC(noJul_mod)
    AIC(Jul_mod)
    AIC(poly_mod)
    AIC(elev_mod_sp)

    pred.y.fall <- predict(mod, newdata = test.data)
    pred.y.fall[pred.y.fall< -0.5] = -0.5
    RMSEP <- sqrt(mean((test.data$y - pred.y.fall)^2))
    SE <- (test.data$y - pred.y.fall)^2
    
    pred.out <- matrix(nrow=length(pred.y.fall), ncol=6)
    pred.out[,1] <- test.data$y
    pred.out[,2] <- pred.y.fall
    pred.out[,3] <- as.character(test.data$SiteName)
    pred.out[,4] <- "fall"
    pred.out[,5] <- q
    pred.out[,6] <- SE
    
    write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD.csv",sep = ",", col.names=F)  
    
    stats.out <- matrix(nrow=1, ncol=10) 
    stats.out[1,1] <- q
    stats.out[1,2] <- "fall"
    stats.out[1,3] <- sum_mod$sigma
    stats.out[1,4] <- sum_mod$adj.r.squared
    stats.out[1,5] <- coeffs[1,1]
    stats.out[1,6] <- coeffs[2,1]
    stats.out[1,7] <- coeffs[3,1]
    stats.out[1,8] <- coeffs[4,1]
    stats.out[1,9] <- coeffs[5,1]
    stats.out[1,10] <- RMSEP

#############All_jk

setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/")

  for (i in 2000:2009)
  {

  spring.data.yr <- spring.data[spring.data$Year == i,]
  fall.data.yr <- fall.data[fall.data$Year == i,]
  
  fall.table <- as.data.frame(table(fall.data.yr$SiteName))  
  SiteID.fall <- as.matrix(fall.table[fall.table$Freq > 10, "Var1"])
  site.num <- nrow(SiteID.fall)
  
  setwd("D:/OneDrive/work/research/Steelhead/Data_analysis/LST/predictions/site_regressions/all_years/jk_all_years")

    j <- 0
    
    while(j <= 100)
    {        
      sites <- sample(SiteID.fall, site.num-1, replace = FALSE, prob=NULL) #change sample size with each data subset
      
      data <- fall.data[fall.data$SiteName %in% sites,]
      test.data <- fall.data[fall.data$SiteName %nin% sites,]
      
      y <- data$y
      x <- data$x
      z <- data$z
      e <- data$e
      mod <- lm(y ~ x + I(x^2) + z + e)
      sum_mod <- summary(mod)
      coeffs <- as.matrix(coefficients(mod))
      
      pred.y.fall <- predict(mod, newdata = test.data)
      pred.y.fall[pred.y.fall< -0.5] = -0.5
      RMSEP <- sqrt(mean((test.data$y - pred.y.fall)^2))
      SE <- (test.data$y - pred.y.fall)^2
      
      pred.out <- matrix(nrow=length(pred.y.fall), ncol=6)
      pred.out[,1] <- test.data$y
      pred.out[,2] <- pred.y.fall
      pred.out[,3] <- as.character(test.data$SiteName)
      pred.out[,4] <- "fall"
      pred.out[,5] <- q
      pred.out[,6] <- SE
      
      write.table (x=pred.out,append=T,row.names=F,file="jk_pred_v_y_JD.csv",sep = ",", col.names=F)  
      
      stats.out <- matrix(nrow=1, ncol=10) 
      stats.out[1,1] <- q
      stats.out[1,2] <- "fall"
      stats.out[1,3] <- sum_mod$sigma
      stats.out[1,4] <- sum_mod$adj.r.squared
      stats.out[1,5] <- coeffs[1,1]
      stats.out[1,6] <- coeffs[2,1]
      stats.out[1,7] <- coeffs[3,1]
      stats.out[1,8] <- coeffs[4,1]
      stats.out[1,9] <- coeffs[5,1]
      stats.out[1,10] <- RMSEP
      
      write.table (x=stats.out, append=T, row.names=F, file="jk_pred_RMSEP_JD.csv", sep=",", col.names=F)
      
      j <- j+1
    }
