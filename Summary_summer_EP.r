#R code to calculate summary metrics for 1July - 30 Sept,  
#27 aug 2015

# first -------------------------------------------------------------------



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



Mean.in <- read.dbf("predt2012_Tuc_8D_Mn_summer.dbf") 
rcas <- unique(Mean.in$RCAID)
SumSumm.out <- data.frame (Mean.in$RCAID)
 

    for (i in 1:length(rcas)) 
        { 
    
        MeanRCA <- Mean.in[i,] #grab days for one RCA 
        DaysAbove12 <- length(which(MeanRCA[1:12]> 12))*8 #finds how many day in the 1July-30 Sept window exceed 13' 
        DaysAbove13 <- length(which(MeanRCA[1:12]> 13))*8
        DaysAbove16 <- length(which(MeanRCA[1:12]> 16))*8
        DaysAbove18 <- length(which(MeanRCA[1:12]> 18))*8
        DaysAbove20 <- length(which(MeanRCA[1:12]> 20))*8
        DaysAbove22 <- length(which(MeanRCA[1:12]> 22))*8
        MaxMean <- max(MeanRCA[1:12])
        SDMean <- sd(MeanRCA[1:12])
        SumSumm.out$CntDays12[i] <- DaysAbove12
        SumSumm.out$CntDays13[i] <- DaysAbove13
        SumSumm.out$CntDays16[i] <- DaysAbove16
        SumSumm.out$CntDays18[i] <- DaysAbove18
        SumSumm.out$CntDays20[i] <- DaysAbove20
        SumSumm.out$CntDays22[i] <- DaysAbove22
        SumSumm.out$MxMn[i] <- MaxMean 
        SumSumm.out$SdMn[i] <- SDMean
        } 
        

SumSumm.out[,8:9] <- round(SumSumm.out[,8:9], digits = 2)
colnames(SumSumm.out)[1] <- "RCAID"
write.dbf(SumSumm.out, file = "Tuc_2012_temp_summary.dbf")

#repeat to calculate summary metrics for 15July - 31Aug,  
#28 aug 2015

# 15Jul-31Aug ------------------------------------------------------------------



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



  Mean.in <- read.dbf("predt2012_Tuc_8D_Mn_summer.dbf") 
  rcas <- unique(Mean.in$RCAID)
  SumSumm.out <- data.frame (Mean.in$RCAID)
  Au.in <- read.csv("Tuc_RCAID_AU.csv")
  Mean.au <- merge(Mean.in, Au.in, by.x = "RCAID", by.y = "rca_id")
  Mean.au$AU_Code <- as.character(Mean.au$AU_Code)
  Aus <- unique(Mean.au$AU_Code)
  
    for (i in 1:length(rcas)) 
      { 
        
        MeanRCA <- Mean.in[i,] #grab days for one RCA 
        DaysAbove12 <- length(which(MeanRCA[2:8]> 12))*8 #finds how many day in the 1July-30 Sept window exceed 13' 
        DaysAbove13 <- length(which(MeanRCA[2:8]> 13))*8
        DaysAbove16 <- length(which(MeanRCA[2:8]> 16))*8
        DaysAbove18 <- length(which(MeanRCA[2:8]> 18))*8
        DaysAbove20 <- length(which(MeanRCA[2:8]> 20))*8
        DaysAbove22 <- length(which(MeanRCA[2:8]> 22))*8
        MaxMean <- max(MeanRCA[1:12])
        SDMean <- sd(MeanRCA[1:12])
        SumSumm.out$CntDays12[i] <- DaysAbove12
        SumSumm.out$CntDays13[i] <- DaysAbove13
        SumSumm.out$CntDays16[i] <- DaysAbove16
        SumSumm.out$CntDays18[i] <- DaysAbove18
        SumSumm.out$CntDays20[i] <- DaysAbove20
        SumSumm.out$CntDays22[i] <- DaysAbove22
        SumSumm.out$MxMn[i] <- MaxMean 
        SumSumm.out$SdMn[i] <- SDMean
      } 

  
  SumSumm.out[,8:9] <- round(SumSumm.out[,8:9], digits = 2)
  colnames(SumSumm.out)[1] <- "RCAID"
  write.dbf(SumSumm.out, file = "Tuc_2012_15Jul_31Aug_summary.dbf")
  write.csv(SumSumm.out, file = "Tuc_2012_15Jul_31Aug_summary.csv", row.names = F)



  MnMn = NULL
  MxMn = NULL
  SdMn = NULL
  AU_Code = NULL

  

    for (i in 1:length(Aus)) 
      { 
        
        MeanAU <- Mean.au[Mean.au$AU_Code == Aus[i],] #grab days for one AU
        MaxMean <- max(MeanAU[,3:9])
        MeanMean <- mean(unlist(MeanAU[3:9]))
        SDMean <- sd(MeanRCA[3:8])
        MnMn = append(MnMn, MeanMean)
        MxMn = append(MxMn, MaxMean)
        SdMn = append(SdMn, SDMean)
        AU_Code = append(AU_Code, Aus[i])
        
      } 

  AuSumm.out = data.frame(MnMn, MxMn, SdMn, AU_Code)
  write.csv(AuSumm.out, file = "Tuc_AU_2012_15Jul_31Aug_summary.csv", row.names = F)


# summary by AU -----------------------------------------------------------



  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "rca_id")
  v12DMn = NULL
  v12Dsd = NULL
  v13DMn = NULL
  v13Dsd = NULL  
  v16DMn = NULL
  v16Dsd = NULL
  v18DMn = NULL
  v18Dsd = NULL
  v20DMn = NULL
  v20Dsd = NULL
  v22DMn = NULL
  v22Dsd = NULL
  AU_Code = NULL


    for (i in 1:length(Aus)) 
      { 
        SumAU <- SumSummAU[SumSummAU$AU_Code == Aus[i],]
        v12Mn = mean(unlist(SumAU$CntDays12))
        v12sd = sd(unlist(SumAU$CntDays12))
        v13Mn = mean(unlist(SumAU$CntDays13))
        v13sd = sd(unlist(SumAU$CntDays13))  
        v16Mn = mean(unlist(SumAU$CntDays16))
        v16sd = sd(unlist(SumAU$CntDays16))
        v18Mn = mean(unlist(SumAU$CntDays18))
        v18sd = sd(unlist(SumAU$CntDays18))
        v20Mn = mean(unlist(SumAU$CntDays20))
        v20sd = sd(unlist(SumAU$CntDays20))
        v22Mn = mean(unlist(SumAU$CntDays22))
        v22sd = sd(unlist(SumAU$CntDays22))
        
        v12DMn = append(v12DMn, v12Mn)
        v12Dsd = append(v12Dsd, v12sd)
        v13DMn = append(v13DMn, v13Mn)
        v13Dsd = append(v13Dsd, v13sd)  
        v16DMn = append(v16DMn, v16Mn)
        v16Dsd = append(v16Dsd, v16sd)
        v18DMn = append(v18DMn, v18Mn)
        v18Dsd = append(v18Dsd, v18sd)
        v20DMn = append(v20DMn, v20Mn)
        v20Dsd = append(v20Dsd, v20sd)
        v22DMn = append(v22DMn, v22Mn)
        v22Dsd = append(v22Dsd, v22sd)
        AU_Code = append(AU_Code, Aus[i])
        
        
      } 
        
    SumAU.out = data.frame(v12DMn, v12Dsd, v13DMn, v13Dsd, v16DMn, v16Dsd, v18DMn, v18Dsd, v20DMn, v20Dsd, v22DMn, v22Dsd, AU_Code)
      
    write.csv(SumAU.out, file = "Tuc_AU_2012_15Jul_31Aug_count_summary.csv", row.names = F)
  
      