############################################################################################################
# These are the model diagnotic tests I use for LST temp models
# Created 14 March 2016, KMcN


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

mainPath <- "D:/OneDrive/work/research/CHaMP/"
basin <- "Ent"
longBasin <- "Entiat"
subDir <- "CHaMP_data"
yrPath1 <- 12
yrPath2 <- 13
yearPath1 <- 2012
yearPath2 <- 2013
setwd(paste0(mainPath, "/", subDir, "/", longBasin))

##########################################################################################
# This includes an example process with a seasonal split to put the diagnostics in context
# The first statement reads in a data table with 5 named columns: x, y, z, e, SiteName.
# These correspond to: LST, Stream temp, Julian Day, Elevation, SiteID
##########################################################################################


  data.in <- read.csv("MyDataSet.csv", stringsAsFactors = FALSE)
  
  split.date <- paste0(yrPath2, "153")
  
  minrow <- which.min(data.in$y)
  data.cool <- data.in[1:minrow,]
  maxrow <- which.max(data.cool$y)
  data.cool <- data.cool[(maxrow+1):dim(data.cool)[1],]
  data.warm <- data.in[minrow:nrow(data.in),]
  data.warm <- data.warm[data.warm$z < split.data,]
  
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

  pred.y[pred.y < -0.5] = -0.5 #this statement sets a floor for water temp predictions (which don't ever go much below 0'C)

  plot(pred.y, y, main = "8-day Min Cooling Leg")
  abline(0,1)

  post_mod <- summary(lm(y ~ pred.y)) #This statement allows me to recalc the RMSE after setting the floor on predictions

  gvmodel <- gvlma(mod)
  summary(gvmodel)
  outlierTest(mod)
  qqPlot(mod, main="QQ Plot Cooling Leg")
  spreadLevelPlot(mod)
  plot(pred.y, mod$residuals, main="Model diagnostics Cooling Leg", xlab="Predicted", ylab="Residuals")

  pressstat_sum <- PRESS(sum_mod, verbose = "FALSE")
  RMSEP <- sqrt(mean((pressstat_sum$residuals)^2))

  library(pls) #this library has to be loaded, run, then unloaded because it masks functions from other libraries
  mod2 <- plsr(y ~ x + I(x^2) + e, validation = "LOO") #Leave-one-out cross validation
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

  mod <- lm(y ~ x + I(x^2) + e)
  sum_mod <- summary(mod)
  coeffs <- as.matrix(coefficients(mod))

  pred.y <- predict(mod)
  pred.y[pred.y < -0.5] = -0.5

  plot(pred.y, y, main = "8-day Min Warming Leg")
  abline(0,1)

  post_mod <- summary(lm(y ~ pred.y))

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

#This section lets me assess estimate quality for the entire dataset

  pred.y <- read.csv(paste0("jk_pred_v_y_Min_", basin, "_", yrPath1, "_", yrPath2, "_sp_fall.csv"), stringsAsFactors = FALSE)
  colnames(pred.y) <- c("Y", "PredY", "JulDay", "Season", "Year")

  plot(pred.y$PredY, pred.y$Y, pch=16, col="blue", main="Min 8-day stream temp 2012-2013", xlab="Predicted", ylab="Observed")
  abline(0,1)
  abline(lm(pred.y$Y~ pred.y$PredY), col="blue")
  fit <- lm(pred.y$Y~ pred.y$PredY)
  plot(fit)


