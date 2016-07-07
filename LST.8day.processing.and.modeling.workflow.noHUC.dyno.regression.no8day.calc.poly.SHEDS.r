############################################################################################################
# This set of R scripts make up the work flow to estimate SHEDS rca temps based on coefficients from standard 
# RCA models
# edited 15 Sept 2015

####################################################################
#This part calculates the zonal mean of the daily LST for each RCA
####################################################################
setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Tucannon/")

out_Preds <- read.dbf("LST11_Tuc_interp.dbf") #use the appropriate read statement

setwd("D:/OneDrive/work/research/CHaMP/CHaMP_data/Tucannon/SHEDS")

weights <- read.csv("SHEDS_rca_area_wgts.csv")
weights[weights$area_wgt<0.01, "area_wgt"] = 0.01

rcas <- unique(unlist(weights$RCA_ID))
rca_zonal <- matrix(nrow = length(rcas), ncol = 46) #change for the number of actual RCAs in the watershed

    l <- 1
    for(i in rcas)  
      {
        pixels <- weights[weights$RCA_ID == i, "GRIDCODE"]
        wgts <- weights[weights$RCA_ID == i, "area_wgt"]
        
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

write.dbf(rca_zonal, file = "LST11_Tuc_SHEDS_RCA.dbf")



########################################################################################################
# This part applies the model coefficients to the LST to generate 8-day temp estimates for 1 July - 30 Sept
########################################################################################################


LST.in <- rca_zonal
colnames(LST.in)<-gsub("X", "", colnames(LST.in))
namesnum <- as.numeric(colnames(LST.in)[1:46])
colnames(LST.in)[1:46] <- namesnum;

elev.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/Tucannon/SHEDS/Tuc_SHEDS_elev.csv")
LST.elev <- merge(LST.in, elev.in, by.x = "RCAID", by.y = "RCA_ID")

coeffs.in <- read.csv("D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/Tuc/2011/All_data_2011_mod_coeffs_Max.csv")
LogPred.out <- LST.elev[,c(25:36,48,1)]
LogPred.out[,1:12] <- 0
rcas <- unique(elev.in$RCA_ID)
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

names.out <- sprintf("Tmn_12_%03d", namesnum)
colnames(LogPred.out)[1:12] <- names.out[1:12]

write.dbf(LogPred.out, file = "D:/OneDrive/work/research/CHaMP/CHaMP_data/Tucannon/SHEDS/predt2011_Tuc_SHEDS_8D_Max_summer.dbf") 

