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



Max.in <- read.dbf("predt2013_Tuc_8D_Max_summer.dbf") 
rcas <- unique(Max.in$RCAID)
SumSumm.out <- data.frame (Max.in$RCAID)
 

    for (i in 1:length(rcas)) 
        { 
    
        MeanRCA <- Max.in[i,] #grab days for one RCA 
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
write.dbf(SumSumm.out, file = "Tuc_2013_temp_summary.dbf")

#repeat to calculate summary metrics for 15July - 31Aug,  
#28 aug 2015

# 15Jul-31Aug ------------------------------------------------------------------



mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/Tuc/"
yearPath <- "2012"


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)

setwd(mainPath)

  Au.in <- read.csv("Tuc_RCAID_AU.csv")
  ST_rcas <- read.csv("STHD_Ext_RCAID_Tuc.csv")
  CH_rcas <- read.csv("CHNK_Ext_RCAID_Tuc.csv")
  
setwd(paste0(mainPath,yearPath))

  Max.in <- read.dbf("predt2012_Tuc_8D_Max_summer.dbf") 
  
  Max.au <- merge(Max.in, Au.in, by.x = "RCAID", by.y = "rca_id")
  Max.au$AU_Code <- as.character(Max.au$AU_Code)
  Aus <- unique(Max.au$AU_Code)

  Max.daily <- Max.in[, c(rep(3:8, each = 8), 13:15)]
  colnames(Max.daily)[1:48] <- 201:248

  rcas <- unique(ST_rcas$RCAID)
  #rcas <- unique(Au.in$rca_id)
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
  #colnames(SumSumm.out)[1] <- "RCAID"
  write.dbf(SumSumm.out, file = "Tuc_2012_21Jul_31Aug_max_summary_STHD.dbf")
  write.csv(SumSumm.out, file = "Tuc_2012_21Jul_31Aug_max_summary_STHD.csv", row.names = F)


# summary by AU -----------------------------------------------------------



  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "rca_id")
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
  write.csv(SumAU.out, file = "Tuc_AU_2012_21Jul_31Aug_ExPcnt_summary_STHD.csv", row.names = F)






# summary by AU including error-----------------------------------------------------------



  SumSummAU <- merge(SumSumm.out, Au.in, by.x = "RCAID", by.y = "rca_id")
  v12DMn = NULL
  v12Dsd = NULL
  v12Dse = NULL
  v13DMn = NULL
  v13Dsd = NULL  
  v13Dse = NULL
  v16DMn = NULL
  v16Dsd = NULL
  v16Dse = NULL
  v18DMn = NULL
  v18Dsd = NULL
  v18Dse = NULL
  v20DMn = NULL
  v20Dsd = NULL
  v20Dse = NULL  
  v22DMn = NULL
  v22Dsd = NULL
  v22Dse = NULL
  MnMxDMn = NULL
  SeMxDMn = NULL
  AU_Code = NULL


    for (i in 1:length(Aus)) 
      { 
        SumAU <- SumSummAU[SumSummAU$AU_Code == Aus[i],]
        v12Mn = mean(unlist(SumAU$PctDays12))
        v12sd = sd(unlist(SumAU$PctDays12))
        v12se = v12sd/sqrt(length(SumAU$PctDays12))
        v13Mn = mean(unlist(SumAU$PctDays13))
        v13sd = sd(unlist(SumAU$PctDays13))  
        v13se = v13sd/sqrt(length(SumAU$PctDays13))
        v16Mn = mean(unlist(SumAU$PctDays16))
        v16sd = sd(unlist(SumAU$PctDays16))
        v16se = v16sd/sqrt(length(SumAU$PctDays16))
        v18Mn = mean(unlist(SumAU$PctDays18))
        v18sd = sd(unlist(SumAU$PctDays18))
        v18se = v18sd/sqrt(length(SumAU$PctDays18))
        v20Mn = mean(unlist(SumAU$PctDays20))
        v20sd = sd(unlist(SumAU$PctDays20))
        v20se = v20sd/sqrt(length(SumAU$PctDays20))
        v22Mn = mean(unlist(SumAU$PctDays22))
        v22sd = sd(unlist(SumAU$PctDays22))
        v22se = v22sd/sqrt(length(SumAU$PctDays22))
        MnMax = mean(unlist(SumAU$MnMx))
        SeMax = MnMax/sqrt(length(SumAU$MnMx))
        
        v12DMn = append(v12DMn, v12Mn)
        v12Dsd = append(v12Dsd, v12sd)
        v12Dse = append(v12Dse, v12se)
        v13DMn = append(v13DMn, v13Mn)
        v13Dsd = append(v13Dsd, v13sd)  
        v13Dse = append(v13Dse, v13se)
        v16DMn = append(v16DMn, v16Mn)
        v16Dsd = append(v16Dsd, v16sd)
        v16Dse = append(v16Dse, v16se)
        v18DMn = append(v18DMn, v18Mn)
        v18Dsd = append(v18Dsd, v18sd)
        v18Dse = append(v18Dse, v18se)
        v20DMn = append(v20DMn, v20Mn)
        v20Dsd = append(v20Dsd, v20sd)
        v20Dse = append(v20Dse, v20se)
        v22DMn = append(v22DMn, v22Mn)
        v22Dsd = append(v22Dsd, v22sd)
        v22Dse = append(v22Dse, v22se)
        MnMxDMn = append(MnMxDMn, MnMax)
        SeMxDMn = append(SeMxDMn, SeMax)
        AU_Code = append(AU_Code, Aus[i])
        
        
      } 
        
    SumAU.out = data.frame(v12DMn, v12Dsd, v12Dse, v13DMn, v13Dsd, v13Dse, v16DMn, v16Dsd, v16Dse, v18DMn, v18Dsd, v18Dse, v20DMn, v20Dsd, v20Dse, v22DMn, v22Dsd, v22Dse, MnMxDMn, SeMxDMn, AU_Code)
    SumAU.out[,1:20] <- round(SumAU.out[,1:20], digits = 2)

    colnames(SumAU.out) <- c("v12DMax", "v12Dsd", "v12Dse", "v13DMax", "v13Dsd", "v13Ds3", "v16DMax", "v16Dsd", "v16Dse", "v18DMax", "v18Dsd", "v18Dse", "v20DMax", "v20Dsd", "v20Dse", "v22DMax", "v22Dsd", "v22Dse", "MnMax", "SeMnMx", "AU_Code")
    write.csv(SumAU.out, file = "Tuc_AU_2013_21Jul_31Aug_ExCnt_summary_STHD.csv", row.names = F)
 






# Degree Days in Exceedence -----------------------------------------------


for (i in 1:length(rcas)) 
{ 
  
  MaxRCA <- Max.in[i,] #grab days for one RCA 
  DaysAbove12 <- as.numeric(MaxRCA[MaxRCA[2:8]> 12][2:8]) #finds days in the 1July-30 Sept window in exceedence of threshold 
  DaysAbove13 <- as.numeric(MaxRCA[MaxRCA[2:8]> 13][2:8])
  DaysAbove16 <- as.numeric(MaxRCA[MaxRCA[2:8]> 16][2:8])
  DaysAbove18 <- as.numeric(MaxRCA[MaxRCA[2:8]> 18][2:8])
  DaysAbove20 <- as.numeric(MaxRCA[MaxRCA[2:8]> 20][2:8])
  DaysAbove22 <- as.numeric(MaxRCA[MaxRCA[2:8]> 22][2:8])
  MaxMax <- max(MaxRCA[1:12])
  SDMax <- sd(MaxRCA[1:12])
  SumSumm.out$DD12[i] <- sum(DaysAbove12 - 12)*8
  SumSumm.out$DD13[i] <- sum(DaysAbove13 - 13)*8
  SumSumm.out$DD16[i] <- sum(DaysAbove16 - 16)*8
  SumSumm.out$DD18[i] <- sum(DaysAbove18 - 18)*8
  SumSumm.out$DD20[i] <- sum(DaysAbove20 - 20)*8
  SumSumm.out$DD22[i] <- sum(DaysAbove22 - 22)*8
  SumSumm.out$MxMx[i] <- MaxMax 
  SumSumm.out$SdMx[i] <- SDMax
} 


SumSumm.out[,8:9] <- round(SumSumm.out[,8:9], digits = 2)
colnames(SumSumm.out)[1] <- "RCAID"
write.dbf(SumSumm.out, file = "Tuc_2013_15Jul_31Aug_max_summary.dbf")
write.csv(SumSumm.out, file = "Tuc_2013_15Jul_31Aug_max_summary.csv", row.names = F)


MnMx = NULL
MxMx = NULL
SdMx = NULL
AU_Code = NULL



  for (i in 1:length(Aus)) 
    { 
      
      MaxAU <- Max.au[Max.au$AU_Code == Aus[i],] #grab days for one AU
      MaxMax <- max(MaxAU[,3:9])
      MeanMax <- mean(unlist(MaxAU[3:9]))
      SDMax <- sd(MaxRCA[3:8])
      MnMx = append(MnMx, MeanMax)
      MxMx = append(MxMx, MaxMax)
      SdMx = append(SdMx, SDMax)
      AU_Code = append(AU_Code, Aus[i])
      
    } 

AuSumm.out = data.frame(MnMx, MxMx, SdMx, AU_Code)
colnames(AuSumm.out) <- c("MnMax", "MaxMax", "SdMax", "AU_Code")
write.csv(AuSumm.out, file = "Tuc_AU_2013_15Jul_31Aug_max_summary.csv", row.names = F)




#_________________________________________________________
# Combining and renaming fields
# Added 8 Oct 2015
#__________________________________________________________

# clean-up ----------------------------------------------------------------


library(timeSeries)
library(lattice)
library(foreign)
library(doBy)
library(qpcR)
library(pls)
library(boot)
library(Hmisc)


  mainPath <- "D:/OneDrive/work/research/CHaMP/CHaMP_data/All_CHaMP/EP_temp/Tuc/"
  yearPath <- "2011"
  
  
  setwd(paste0(mainPath,yearPath))
  
  Max.in <- read.dbf("Tuc_2011_21Jul_31Aug_max_summary.dbf") 
  Max.out <- Max.in
  
  colnames(Max.out) <- c("RCAID", "Pct12_2011","Pct13_2011", "Pct16_2011", "Pct18_2011", "Pct20_2011", "Pct22_2011",  "MxMx_2011", "sdMn_2011", "MnMx_2011")
  
  write.dbf(Max.out, file = "Tuc_2011_21Jul_31Aug_max_summary.dbf")










