############################################################################################################
# This set of R scripts generates 8D temp estimates for July-Sept only
# created 13 July 2015 
# Tucannon
          
mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/Tuc"
setwd(mainPath)

library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

#################################################################
# 2012
##############
# copied the NoNA.xyz from the SA analysis directory to use here

NoNA.xyz <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/Sensitivity_analysis/Tucannon_2012/NoNA_xyz.csv")

coeffs_out <- data.frame(Int=numeric(5), bLST=numeric(5), bLST2=numeric(5), bJul=numeric(1), bElev=numeric(5))
metrics_out <- data.frame(PRESS=numeric(5), p2=numeric(5), r2=numeric(5), RMSEP=numeric(5), RMSE=numeric(5), AIC=numeric(5))
rownames(metrics_out) <- c("noJul", "Jul", "Poly", "Elev", "All")
rownames(coeffs_out) <- c("noJul", "Jul", "Poly", "Elev", "All")

    summer <- subset(NoNA.xyz, JulDay >181 & JulDay < 274)  
   

     y <- summer$Temp8DayMax
     x <- summer$LST
     z <- summer$JulDay
     e <- summer$Elev
      plot(x, y)
      
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

        pred.y.spring[pred.y.spring < -0.5] = -0.5
        pred.out <- matrix(nrow=length(y), ncol=4)
        colnames(pred.out) <- c("Y_max", "PredY_max", "SiteName", "Year")
        pred.out[,1] <- y
        pred.out[,2] <- pred.y.spring
        pred.out[,3] <- as.character(summer$SiteName)
        pred.out[,4] <- "Tuc"
        plot(pred.out[,1], pred.out[,2])
        
     
      write.table (x=pred.out,append=F,row.names=F,file="2012_pred_v_y_max_Tuc.csv",sep = ",", col.names=T)  
      
     
     
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


write.table(x=coeffs_out, append=F,row.names=F, file = "All_data_2012_mod_coeffs_Max.csv", sep = ",", col.names=T)

write.table(x=metrics_out, append=F,row.names=F, file = "All_data_2012_mod_metrics_Max.csv", sep = ",", col.names=T)






  


##############################################
# Full model data jk N-1 for stats
# Takes the full data sets and jk out 1 site for validation for 100 reps 

# jackknife ---------------------------------------------------------------




pathname <- "D:\\OneDrive\\work\\research\\CHaMP\\CHaMP_data\\All_CHaMP\\EP_temp\\Tuc\\jk_for_val\\"
setwd(pathname)

pred.out <- matrix(nrow=0, ncol = 5)
colnames(pred.out) <- c("Y_sp", "PredY_sp", "SiteName", "Year", "Basin")
write.table (x=pred.out,append=F,row.names=F,file="2012_jk_pred_v_y_Tuc.csv",sep = ",", col.names=T)

j <- 0

  while(j <= 100)
    {        
      sites <- sample(SiteID.sp, 15, replace = FALSE, prob=NULL) #change sample size with each data subset
      
      write.table(x=sites, append=T,row.names=F, file = "sites_jk_n-1.csv", sep = ",", col.names=F)
      
      data <- summer[summer$SiteName %in% sites,]
      test.data <- summer[summer$SiteName %nin% sites,]
      
      y <- data$y
      x <- data$x
      z <- data$z
      e <- data$e
      mod <- lm(y ~ x + I(x^2) + z + e)
      
      pred.y.spring <- predict(mod, newdata = test.data)
      
      pred.y.spring[pred.y.spring<0.5] = 0
      pred.out <- matrix(nrow=length(pred.y.spring), ncol=5)
      pred.out[,1] <- test.data$y
      pred.out[,2] <- pred.y.spring
      pred.out[,3] <- as.character(test.data$SiteName)
      pred.out[,4] <- "2012"
      pred.out[,5] <- "Tuc"
      #plot(pred.out[,1], pred.out[,2])
      write.table (x=pred.out,append=T,row.names=F,file="2012_jk_pred_v_y_Tuc.csv",sep = ",", col.names=F)  
      
      j <- j+1
    }




##############################

########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day MAX temp estimates for 1 July - 30 Sept
########################################################################################################


setwd(mainPath)
LST.in <- read.dbf("D:/Dropbox/work/research/CHaMP/CHaMP_data/Tucannon/LST12_Tuc_RCA.dbf")
colnames(LST.in)<-gsub("X", "", colnames(LST.in))
namesnum <- as.numeric(colnames(LST.in)[1:46])
colnames(LST.in)[1:46] <- namesnum;

elev.in <- read.csv("D:/Dropbox/work/research/CHaMP/CHaMP_data/Tucannon/Tuc_PtID_elev.csv")
LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "GRIDCODE")

coeffs.in <- read.csv("All_data_2012_mod_coeffs_Max.csv")
LogPred.out <- LST.elev[,c(25:36,48,1)]
LogPred.out[,1:12] <- 0
rcas <- unique(elev.in$GRIDCODE)
LST.summ <- LST.elev[,c(25:36, 48,1)]

#############################
#basin-wide model application
#############################

  
  for (i in 1:length(rcas))  
    {
      x <- unlist(LST.summ[i,])
      
      
      j <- 185
      
      for (l in 1:12)
      {x[l] <- x[l] * coeffs.in$bLST[5] + j * coeffs.in$bJul[5] +  x[13] * coeffs.in$bElev[5] + x[l]^2 * coeffs.in$bLST2[5] + coeffs.in$Int[5]
       j <- j + 8}
     
      LogPred.out[i,1:12] <- x [1:12] 
    }

  
LogPred.out <- as.data.frame(LogPred.out)
LogPred.out[LogPred.out< -0.5] = -0.5
LogPred.out$Basin_RCA <- paste0("Tucannon_", LogPred.out$RCAID)
namesnum <- as.numeric(colnames(LogPred.out[1:12]))

names.out <- sprintf("Tmx_12_%03d", namesnum)
colnames(LogPred.out)[1:12] <- names.out[1:12]

write.dbf(LogPred.out, file = "predt2012_Tuc_8D_Max_summer.dbf") 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  