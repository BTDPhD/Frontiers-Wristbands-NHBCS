##############################

#Title: Assessment of Multipollutant Exposures during Pregnancy using Silicone Wristbands
#Authors: Brett T Doherty, John L Pearce, Kim A Anderson, Margaret R Karagas, Megan E Romano

##############################

#This file includes the code to prepare and analyze the data pertaining to the 24 gestational week wristbands

##############################

#Contents

#AAA - Libraries
#BBB - Data
#CCC - Play Area

#DDD - Data Manipulations
  #DDD1 - Create "Chem" vars
  #DDD2 - Create "Chem_det" vars
  #DDD3 - Create total detect count variable
  #DDD4_12W - Create category-specific detect count variables
  #DDD4_24W - Create category-specific detect count variables
  #DDD6 - Covariates
  #DDD7 - Add 12W standardized variables
  #DDD71 - Add 24W standardized variables
  #DDD11 - Add 12W batch-corrected centered-scaled
  #DDD111 - Add 24W batch-corrected centered-scaled

#EEE - Table 1 - Population Characteristics
#FFF - Table 2 - Summary Statistics - Wristbands
#III - Wristband-Wristband Relations - ICCs

##############################

#AAA - Libraries

library(tidyverse)

##############################

#BBB - Data

#Set directory
#XXX

#Load data
#XXX

#Set data to be cleaned
toClean <- clean

##############################

#DDD - Data Manipulations

#DDD1 - Create "Chem" vars

#Create dataset that contains only the 12W CASN variables
CASNvarsold <- select(toClean, starts_with("M_"))
CASNvarsold12W <- select(CASNvarsold, ends_with("12W"))
CASNvarsold24W <- select(CASNvarsold, ends_with("24W"))

#Sort dataset to ensure numeric/alphabetical order
CASNvarsold12W<-CASNvarsold12W[,order(colnames(CASNvarsold12W),decreasing=FALSE)]
CASNvarsold24W<-CASNvarsold24W[,order(colnames(CASNvarsold24W),decreasing=FALSE)]

#Create column matrix of names of old CASN variables
CASNvarsold12WNames <- as.matrix(names(CASNvarsold12W))
CASNvarsold24WNames <- as.matrix(names(CASNvarsold24W))

#Create column matrix of names for new CASN variables
CASNvarsnewnames12W <- matrix(NA, nrow = length(CASNvarsold12WNames), ncol = 1)
for(i in 1:length(CASNvarsold12WNames)){CASNvarsnewnames12W[i] <- paste("chem", i,"12W", sep = "_")}

CASNvarsnewnames24W <- matrix(NA, nrow = length(CASNvarsold24WNames), ncol = 1)
for(i in 1:length(CASNvarsold24WNames)){CASNvarsnewnames24W[i] <- paste("chem", i,"24W", sep = "_")}

#Create matrix that includes both, to be used as an index
chemBook <- matrix(NA, nrow = length(CASNvarsold12WNames), ncol = 4)
chemBook[,1] <- CASNvarsold12WNames
chemBook[,2] <- CASNvarsnewnames12W
chemBook[,3] <- CASNvarsold24WNames
chemBook[,4] <- CASNvarsnewnames24W

#XXX

#Apply new names to old dataset
CASNvarsnew12W <- CASNvarsold12W
for(i in 1:ncol(CASNvarsnew12W)){colnames(CASNvarsnew12W)[i] <- CASNvarsnewnames12W[i]}

CASNvarsnew24W <- CASNvarsold24W
for(i in 1:ncol(CASNvarsnew24W)){colnames(CASNvarsnew24W)[i] <- CASNvarsnewnames24W[i]}

#Add these variables to the larger dataset
const <- cbind(toClean, CASNvarsnew12W, CASNvarsnew24W)

#Write to wrist
#XXX

##############################

#DDD2 - Create Chem Det Vars

#Take new var dataset
CASNvarsnew12WDet <- CASNvarsnew12W
CASNvarsnew24WDet <- CASNvarsnew24W

#Apply function to each column that creates a new variable: 1 if  chemX != 0, 0 if chemX = 0
CASNvarsnew12WDet[, paste0(names(CASNvarsnew12WDet),"_det")] <- lapply(CASNvarsnew12WDet, function(x) as.numeric(x!=0))
CASNvarsnew24WDet[, paste0(names(CASNvarsnew24WDet),"_det")] <- lapply(CASNvarsnew24WDet, function(x) as.numeric(x!=0))

#Create dataset with these variables only
detVars12W <- select(CASNvarsnew12WDet, ends_with("_det"))
detVars24W <- select(CASNvarsnew24WDet, ends_with("_det"))

#Add these variables to the largerdataset
const <- cbind(const, detVars12W, detVars24W)

#Write to wrist
#XXX

#DDD3 - Create total detect count variable

#Create new variable, Total, which is the sum of the det vars for each participant
detVars12W$Total_12W = rowSums(detVars12W)
detVars24W$Total_24W = rowSums(detVars24W)

#Isolate this variable in preparation to append to larger dataset
totalVar12W <- select(detVars12W, Total_12W)
totalVar24W <- select(detVars24W, Total_24W)

#Add variable to the larger dataset
const <- cbind(const, totalVar12W, totalVar24W)

#Write to wrist
#XXX

##############################

#DDD4_12W - Create category-specific detect count variables

#XXX
#XXX
#XXX

niceClass <- class

#Make sets of each of the variables by their classifications
PAH <- niceClass[which(niceClass$PAH == 1),]
PAH_names <- as.matrix(PAH$chem_12W)

FR <- niceClass[which(niceClass$FR == 1),]
FR_names <- as.matrix(FR$chem_12W)

Phar <- niceClass[which(niceClass$Phar == 1),]
Phar_names <- as.matrix(Phar$chem_12W)

PC <- niceClass[which(niceClass$PC == 1),]
PC_names <- as.matrix(PC$chem_12W)

CiC <- niceClass[which(niceClass$CiC == 1),]
CiC_names <- as.matrix(CiC$chem_12W)

CP <- niceClass[which(niceClass$CP == 1),]
CP_names <- as.matrix(CP$chem_12W)

Pest <- niceClass[which(niceClass$Pest == 1),]
Pest_names <- as.matrix(Pest$chem_12W)

#Sum _det variables over each class;

#CiC
CiC_names_all <- matrix(NA, nrow = (length(CiC_names)), ncol = 1)
for(i in 1:length(CiC_names)){CiC_names_all[i,] <- paste(CiC_names[i],"_det", sep = "")}
CiCVarsDet <- const[,CiC_names_all]
CiCVarsDet$TotalCiC_12W = rowSums(CiCVarsDet)
TotalCiC_12W <- select(CiCVarsDet, TotalCiC_12W)

#CP
CP_names_all <- matrix(NA, nrow = (length(CP_names)), ncol = 1)
for(i in 1:length(CP_names)){CP_names_all[i,] <- paste(CP_names[i],"_det", sep = "")}
CPVarsDet  <- const[,CP_names_all]
CPVarsDet$TotalCP_12W = rowSums(CPVarsDet)
TotalCP_12W <- select(CPVarsDet, TotalCP_12W)

#FR
FR_names_all <- matrix(NA, nrow = (length(FR_names)), ncol = 1)
for(i in 1:length(FR_names)){FR_names_all[i,] <- paste(FR_names[i],"_det", sep = "")}
FRVarsDet  <- const[,FR_names_all]
FRVarsDet$TotalFR_12W = rowSums(FRVarsDet)
TotalFR_12W <- select(FRVarsDet, TotalFR_12W)

#PAH
PAH_names_all <- matrix(NA, nrow = (length(PAH_names)), ncol = 1)
for(i in 1:length(PAH_names)){PAH_names_all[i,] <- paste(PAH_names[i],"_det", sep = "")}
PAHVarsDet  <- const[,PAH_names_all]
PAHVarsDet$TotalPAH_12W = rowSums(PAHVarsDet)
TotalPAH_12W <- select(PAHVarsDet, TotalPAH_12W)

#PC
PC_names_all <- matrix(NA, nrow = (length(PC_names)), ncol = 1)
for(i in 1:length(PC_names)){PC_names_all[i,] <- paste(PC_names[i],"_det", sep = "")}
PCVarsDet  <- const[,PC_names_all]
PCVarsDet$TotalPC_12W = rowSums(PCVarsDet)
TotalPC_12W <- select(PCVarsDet, TotalPC_12W)

#Pest
Pest_names_all <- matrix(NA, nrow = (length(Pest_names)), ncol = 1)
for(i in 1:length(Pest_names)){Pest_names_all[i,] <- paste(Pest_names[i],"_det", sep = "")}
PestVarsDet <- const[,Pest_names_all]
PestVarsDet$TotalPest_12W = rowSums(PestVarsDet)
TotalPest_12W <- select(PestVarsDet, TotalPest_12W)

#Pharm
Phar_names_all <- matrix(NA, nrow = (length(Phar_names)), ncol = 1)
for(i in 1:length(Phar_names)){Phar_names_all[i,] <- paste(Phar_names[i],"_det", sep = "")}
PharVarsDet <- const[,Phar_names_all]
PharVarsDet$TotalPhar_12W = rowSums(PharVarsDet)
TotalPhar_12W <- select(PharVarsDet, TotalPhar_12W)

const <- cbind(const, TotalCiC_12W, TotalCP_12W, TotalFR_12W, TotalPAH_12W, TotalPC_12W, TotalPest_12W, TotalPhar_12W)

#XXX

##############################

#DDD4_24W - Create category-specific detect count variables

#XXX
#XXX
#XXX

niceClass <- class

#Make sets of each of the variables by their classifications
PAH <- niceClass[which(niceClass$PAH == 1),]
PAH_names <- as.matrix(PAH$chem_24W)

FR <- niceClass[which(niceClass$FR == 1),]
FR_names <- as.matrix(FR$chem_24W)

Phar <- niceClass[which(niceClass$Phar == 1),]
Phar_names <- as.matrix(Phar$chem_24W)

PC <- niceClass[which(niceClass$PC == 1),]
PC_names <- as.matrix(PC$chem_24W)

CiC <- niceClass[which(niceClass$CiC == 1),]
CiC_names <- as.matrix(CiC$chem_24W)

CP <- niceClass[which(niceClass$CP == 1),]
CP_names <- as.matrix(CP$chem_24W)

Pest <- niceClass[which(niceClass$Pest == 1),]
Pest_names <- as.matrix(Pest$chem_24W)

#Sum _det variables over each class;

#CiC
CiC_names_all <- matrix(NA, nrow = (length(CiC_names)), ncol = 1)
for(i in 1:length(CiC_names)){CiC_names_all[i,] <- paste(CiC_names[i],"_det", sep = "")}
CiCVarsDet <- const[,CiC_names_all]
CiCVarsDet$TotalCiC_24W = rowSums(CiCVarsDet)
TotalCiC_24W <- select(CiCVarsDet, TotalCiC_24W)

#CP
CP_names_all <- matrix(NA, nrow = (length(CP_names)), ncol = 1)
for(i in 1:length(CP_names)){CP_names_all[i,] <- paste(CP_names[i],"_det", sep = "")}
CPVarsDet <- const[,CP_names_all]
CPVarsDet$TotalCP_24W = rowSums(CPVarsDet)
TotalCP_24W <- select(CPVarsDet, TotalCP_24W)

#FR
FR_names_all <- matrix(NA, nrow = (length(FR_names)), ncol = 1)
for(i in 1:length(FR_names)){FR_names_all[i,] <- paste(FR_names[i],"_det", sep = "")}
FRVarsDet <- const[,FR_names_all]
FRVarsDet$TotalFR_24W = rowSums(FRVarsDet)
TotalFR_24W <- select(FRVarsDet, TotalFR_24W)

#PAH
PAH_names_all <- matrix(NA, nrow = (length(PAH_names)), ncol = 1)
for(i in 1:length(PAH_names)){PAH_names_all[i,] <- paste(PAH_names[i],"_det", sep = "")}
PAHVarsDet <- const[,PAH_names_all]
PAHVarsDet$TotalPAH_24W = rowSums(PAHVarsDet)
TotalPAH_24W <- select(PAHVarsDet, TotalPAH_24W)

#PC
PC_names_all <- matrix(NA, nrow = (length(PC_names)), ncol = 1)
for(i in 1:length(PC_names)){PC_names_all[i,] <- paste(PC_names[i],"_det", sep = "")}
PCVarsDet <- const[,PC_names_all]
PCVarsDet$TotalPC_24W = rowSums(PCVarsDet)
TotalPC_24W <- select(PCVarsDet, TotalPC_24W)

#Pest
Pest_names_all <- matrix(NA, nrow = (length(Pest_names)), ncol = 1)
for(i in 1:length(Pest_names)){Pest_names_all[i,] <- paste(Pest_names[i],"_det", sep = "")}
PestVarsDet <- const[,Pest_names_all]
PestVarsDet$TotalPest_24W = rowSums(PestVarsDet)
TotalPest_24W <- select(PestVarsDet, TotalPest_24W)

#Phar
Phar_names_all <- matrix(NA, nrow = (length(Phar_names)), ncol = 1)
for(i in 1:length(Phar_names)){Phar_names_all[i,] <- paste(Phar_names[i],"_det", sep = "")}
PharVarsDet <- const[,Phar_names_all]
PharVarsDet$TotalPhar_24W = rowSums(PharVarsDet)
TotalPhar_24W <- select(PharVarsDet, TotalPhar_24W)

const <- cbind(const, TotalCiC_24W, TotalCP_24W, TotalFR_24W, TotalPAH_24W, TotalPC_24W, TotalPest_24W, TotalPhar_24W)

#XXX

##############################

#DDD6 - Covariates

#Non-Hispanic White
const$mth_WNH <- ifelse(const$mth_race == "e. White" & const$mth_hisp == 0, 1, 0)
const$mth_WNH_flip <- ifelse(const$mth_race == "e. White" & const$mth_hisp == 0, 0, 1)

#2-level PNQRelatStat
const$PNQRelatStat2Lvl <- ifelse(const$PNQRelatStat == 1, 1, 0)
const$PNQRelatStat2Lvl_flip <- ifelse(const$PNQRelatStat == 1, 0, 1)

#3-level PNQSchoolLvl
const$PNQSchoolLvl3Lvl <- ifelse(const$PNQSchoolLvl == 1 | const$PNQSchoolLvl == 2 | const$PNQSchoolLvl == 3, 0, const$PNQSchoolLvl)

#2-level parity
const$parity2Lvl <- ifelse(const$parity == 0, 0,
                           ifelse(const$parity >= 1, 1, NA))

const$parity2Lvl_flip <- ifelse(const$parity == 0, 1,
                           ifelse(const$parity >= 1, 0, NA))

#3-Level BMI
const$BMI3Lvl <- NA
const$BMI3Lvl <- ifelse(const$maternalbmi == "NA", NA, 
                        ifelse(const$maternalbmi < 18.5, 0, 
                               ifelse(const$maternalbmi < 24.9, 1,  
                                      ifelse(const$maternalbmi < 29.9, 2, 3)))) 

#Two-Level Smoking
const$smokegood <- ifelse(const$evercigpreg == 1, 1, 
                          ifelse(const$shsmokepreg == 1, 1, 
                                 ifelse(const$shsmokepreg == 0 & const$evercigpreg == 0, 0, NA)))
const$smokegood_flip <- ifelse(const$evercigpreg == 1, 0, 
                          ifelse(const$shsmokepreg == 1, 0, 
                                 ifelse(const$shsmokepreg == 0 & const$evercigpreg == 0, 1, NA)))

#Season 12W
const$seasonwear12W <- ifelse(is.na(const$month_start_wear_12W), NA,
                              ifelse(const$month_start_wear_12W %in% c(12,1,2), 1, 
                                     ifelse(const$month_start_wear_12W %in% c(3,4,5), 2,
                                            ifelse(const$month_start_wear_12W %in% c(6,7,8), 3,
                                                   ifelse(const$month_start_wear_12W %in% c(9, 10, 11), 4, 999)))))
const$seasonwear12WWint <- ifelse(is.na(const$seasonwear12W), NA, 
                                  ifelse(const$seasonwear12W == 1, 1, 0))

const$seasonwear12WSpri <- ifelse(is.na(const$seasonwear12W), NA, 
                                  ifelse(const$seasonwear12W == 2, 1, 0))

const$seasonwear12WSumm <- ifelse(is.na(const$seasonwear12W), NA, 
                                  ifelse(const$seasonwear12W == 3, 1, 0))

const$seasonwear12WFall <- ifelse(is.na(const$seasonwear12W), NA, 
                                  ifelse(const$seasonwear12W == 4, 1, 0))

#Season 24W
const$seasonwear24W <- ifelse(is.na(const$month_start_wear_24W), NA,
                              ifelse(const$month_start_wear_24W %in% c(12,1,2), 1, 
                                     ifelse(const$month_start_wear_24W %in% c(3,4,5), 2,
                                            ifelse(const$month_start_wear_24W %in% c(6,7,8), 3,
                                                   ifelse(const$month_start_wear_24W %in% c(9, 10, 11), 4, 999)))))

#WB Size 2-level
const$WBSize2Lvl12W <- ifelse(const$Wristband_Size_12W == "regular", 1,
                              ifelse(const$Wristband_Size_12W == "small", 0, NA))
const$WBSize2Lvl12W_flip <- ifelse(const$Wristband_Size_12W == "regular", 0,
                              ifelse(const$Wristband_Size_12W == "small", 1, NA))

#Flag for 12W
const$OIF_12W_Flag <- ifelse(is.na(const$M_1019_39765_80_5_12W),0,1)

#Flag for 24W
const$OIF_24W_Flag <- ifelse(is.na(const$M_1019_39765_80_5_24W),0,1)

#Flag for both 12W and 24W
const$OIF_12Wand24W_Flag <- ifelse(const$OIF_12W_Flag == 1 & const$OIF_24W_Flag,1,0)

#Flag for either 12W or 24W
const$OIF_12Wor24W_Flag <- ifelse(const$OIF_12W_Flag == 1 | const$OIF_24W_Flag,1,0)

#Gestage Start in Weeks
const$gestage_wbstart_12w_weeks <- const$gestage_wbstart_12w / 7
const$gestage_wbstart_24w_weeks <- const$gestage_wbstart_24w / 7

#Write to CSV
#XXX

##############################

#DDD7 - Add 12W standardized variables

#Step 1 - Create Matrix of batch-chem specific medians

chems <- select(const, starts_with("chem"))
chems <- select(chems, -contains("det"))
chems <- select(chems, -contains("log"))
chems <- select(chems, contains("12W"))

batch <- as.factor(const$Lab_Submission_Batch_12W)

chemsBatch <- cbind(batch, chems)

capMat <- matrix(nrow = ncol(chems), ncol = 5)

for(i in 2:ncol(chemsBatch)){
  
  capMat[i-1,] <- with(chemsBatch, tapply(chemsBatch[,i], chemsBatch[,1], median))
  
}

batchChemMeds <- capMat

#Step 2. Subtract the correct batch-chem median from each value.

table(const$Lab_Submission_Batch_12W)

chems <- select(const, starts_with("chem"))
chems <- select(chems, -contains("det"))
chems <- select(chems, -contains("log"))
chems <- select(chems, contains("12W"))

batch <- as.factor(const$Lab_Submission_Batch_12W)

chemsBatchOrig <- cbind(batch, chems)

chemsBatch <- cbind(chems,batch)

for(i in 1:(ncol(chemsBatch) - 1)){
  for(j in 1:nrow(chemsBatch)){
    
    chemsBatch[j,i] <- chemsBatch[j,i] - batchChemMeds[i, chemsBatch[j,159]]
    
  }
}

#remove the batch var
chemsBatch <- chemsBatch[,1:158]

#Rename the variables
colnames(chemsBatch) <- paste(colnames(chemsBatch), "_std", sep = "")

#Add to toCleanData
const <- cbind(const, chemsBatch)

#Write to CSV
#XXX

##############################

#DDD71 - Add 24W standardized variables

#Plan:
#Create a matrix of all the batch-chemical medians
#Somehow cycle through the chemical matrix, deducting the appropriate value from each

#Step 1 - Create Matrix of batch-chem specific medians

table(const$Lab_Submission_Batch_24W)

chems <- select(const, starts_with("chem"))
chems <- select(chems, -contains("det"))
chems <- select(chems, -contains("log"))
chems <- select(chems, contains("24W"))

batch <- as.factor(const$Lab_Submission_Batch_24W)

chemsBatch <- cbind(batch, chems)

capMat <- matrix(nrow = ncol(chems), ncol = 5)

for(i in 2:ncol(chemsBatch)){
  
  capMat[i-1,] <- with(chemsBatch, tapply(chemsBatch[,i], chemsBatch[,1], median))
  
}

batchChemMeds <- capMat

#Step 2. Subtract the correct batch-chem median from each value.

table(const$Lab_Submission_Batch_24W)

chems <- select(const, starts_with("chem"))
chems <- select(chems, -contains("det"))
chems <- select(chems, -contains("log"))
chems <- select(chems, contains("24W"))

batch <- as.factor(const$Lab_Submission_Batch_24W)

chemsBatchOrig <- cbind(batch, chems)

chemsBatch <- cbind(chems,batch)

for(i in 1:(ncol(chemsBatch) - 1)){
  for(j in 1:nrow(chemsBatch)){
    
    chemsBatch[j,i] <- chemsBatch[j,i] - batchChemMeds[i, chemsBatch[j,159]]
    
  }
}

#remove the batch var
chemsBatch <- chemsBatch[,1:158]

#Rename the variables
colnames(chemsBatch) <- paste(colnames(chemsBatch), "_std", sep = "")

#Add to toCleanData
const <- cbind(const, chemsBatch)

#Write to CSV
#XXX

##############################

#DDD11 - Add 12W batch-corrected centered-scaled

batchStan <- dplyr::select(const, ends_with("std"))
batchStan <- dplyr::select(batchStan, -contains("log2"))
batchStan <- dplyr::select(batchStan, contains("12W"))

batchStanStan <- as.data.frame(scale(batchStan, center = T, scale = T))

#Rename the variables
colnames(batchStanStan) <- paste(colnames(batchStanStan), "_cs", sep = "")

#Add to toCleanData
const <- cbind(const, batchStanStan)

#Write to CSV
#XXX

##############################

#DDD111 - Add 24W batch-corrected centered-scaled

batchStan <- dplyr::select(const, ends_with("std"))
batchStan <- dplyr::select(batchStan, -contains("log2"))
batchStan <- dplyr::select(batchStan, contains("24W"))

batchStanStan <- as.data.frame(scale(batchStan, center = T, scale = T))

#Rename the variables
colnames(batchStanStan) <- paste(colnames(batchStanStan), "_cs", sep = "")

#Add to toCleanData
const <- cbind(const, batchStanStan)

#Write to CSV
#XXX

##############################
##############################
##############################
##############################
##############################

#EEE - Table 1

#Set directory
#XXX

#Load cleaned data
#XXX

#Check samples

table(const$OIF_12W_Flag)
table(const$OIF_24W_Flag)
table(const$OIF_12Wand24W_Flag)
table(const$OIF_12Wor24W_Flag)

#24W Only
OIF24 <- const[which(const$OIF_24W_Flag == 1),]

library(janitor)

summary(OIF24$enrollment_age)
summary(OIF24$maternalbmi)
tabyl(OIF24$PNQSchoolLvl3Lvl)
tabyl(OIF24$PNQRelatStat2Lvl)
tabyl(OIF24$mth_WNH)
tabyl(OIF24$parity2Lvl)
tabyl(OIF24$smokegood)

summary(OIF24$gestage_wbstart_24w_weeks)
summary(OIF24$Days_worn_24W)
tabyl(OIF24$Wristband_Size_24W)
tabyl(OIF24$seasonwear24W)

##############################

#FFF - Summary Statistics - Wristbands

#Set directory
#XXX

#Read cleaned data
#XXX

#Read cleaned classes file
#XXX

#Summary detection for wristbands - 24 week (i.e., % detected for each chemical)

chems <- select(const, contains("chem"))
chems24W <- select(chems, contains("24W"))
chems24Wdet <- select(chems24W, contains("det"))
chems24Wdet <- chems24Wdet[complete.cases(chems24Wdet),]

DET<-function (x){length(which(x == 1))}
DETs<-round(apply(chems24Wdet,2,FUN=DET),2)
DETs.per <- round(100*(DETs/dim(chems24Wdet)[1]),2)
write.csv(DETs.per, "C:\\R_Out\\20200120.csv")

#Basic descriptives for wristbands - 24 week (i.e., n and percentiles for each chemical)

chems <- select(const, contains("chem"))
chems24W <- select(chems, contains("24W"))
chems24Wconc <- select(chems24W, -contains("det"))
chems24Wconc <- select(chems24Wconc, -contains("log"))
chems24Wconc <- select(chems24Wconc, -contains("std"))
chems24Wconc <- chems24Wconc[complete.cases(chems24Wconc),]

p50 <- round(apply(chems24Wconc,2,FUN='quantile',probs = 0.50, na.rm=TRUE),0)
p25 <- round(apply(chems24Wconc,2,FUN='quantile',probs = 0.25, na.rm=TRUE),0)
p75 <- round(apply(chems24Wconc,2,FUN='quantile',probs = 0.75, na.rm=TRUE),0)
p95 <- round(apply(chems24Wconc,2,FUN='quantile',probs = 0.95, na.rm=TRUE),0)
max <- round(apply(chems24Wconc,2,FUN='max', na.rm=TRUE),0)

chems24Wstats <- cbind(p50, p25, p75, p95, max)

#XXX

########################################################

#GGG - Chemicals Summary

#Set directory
#XXX

#Read cleaned data
#XXX

#Restrict to participants at 24 weeks
OIF24 <- const[which(const$OIF_24W_Flag == 1),]

#24 Week

capMat <- matrix(nrow = 8, ncol = 6)

a <- summary(OIF24$Total_24W)
capMat[1,1] <- "Total"
capMat[1,2] <- a[3]
capMat[1,3] <- a[2]
capMat[1,4] <- a[5]
capMat[1,5] <- a[1]
capMat[1,6] <- a[6]

a <- summary(OIF24$TotalCiC_24W)
capMat[2,1] <- "CiC"
capMat[2,2] <- a[3]
capMat[2,3] <- a[2]
capMat[2,4] <- a[5]
capMat[2,5] <- a[1]
capMat[2,6] <- a[6]

a <- summary(OIF24$TotalPC_24W)
capMat[3,1] <- "PC"
capMat[3,2] <- a[3]
capMat[3,3] <- a[2]
capMat[3,4] <- a[5]
capMat[3,5] <- a[1]
capMat[3,6] <- a[6]

a <- summary(OIF24$TotalPest_24W)
capMat[4,1] <- "Pest"
capMat[4,2] <- a[3]
capMat[4,3] <- a[2]
capMat[4,4] <- a[5]
capMat[4,5] <- a[1]
capMat[4,6] <- a[6]

a <- summary(OIF24$TotalFR_24W)
capMat[5,1] <- "FR"
capMat[5,2] <- a[3]
capMat[5,3] <- a[2]
capMat[5,4] <- a[5]
capMat[5,5] <- a[1]
capMat[5,6] <- a[6]

a <- summary(OIF24$TotalPAH_24W)
capMat[6,1] <- "PAH"
capMat[6,2] <- a[3]
capMat[6,3] <- a[2]
capMat[6,4] <- a[5]
capMat[6,5] <- a[1]
capMat[6,6] <- a[6]

a <- summary(OIF24$TotalCP_24W)
capMat[7,1] <- "CP"
capMat[7,2] <- a[3]
capMat[7,3] <- a[2]
capMat[7,4] <- a[5]
capMat[7,5] <- a[1]
capMat[7,6] <- a[6]

a <- summary(OIF24$TotalPhar_24W)
capMat[8,1] <- "Pharm"
capMat[8,2] <- a[3]
capMat[8,3] <- a[2]
capMat[8,4] <- a[5]
capMat[8,5] <- a[1]
capMat[8,6] <- a[6]

#XXX

############################################

#III - Wristband-Wristband Relations - ICCs

library(dplyr)

#Set directory
#XXX

#Read cleaned data
#XXX

#Restrict to participants with both 12 and 24 week data
OIF12and24 <- const[which(const$OIF_12Wand24W_Flag == 1),]

#Select chemicals that have >60% detection at either time point

chems12 <- dplyr::select(OIF12and24,c(
  "chem_137_12W_std_cs",
  "chem_99_12W_std_cs",
  "chem_66_12W_std_cs",
  "chem_17_12W_std_cs",
  "chem_114_12W_std_cs",
  "chem_111_12W_std_cs",
  "chem_67_12W_std_cs",
  "chem_29_12W_std_cs",
  "chem_40_12W_std_cs",
  "chem_100_12W_std_cs",
  "chem_63_12W_std_cs",
  "chem_80_12W_std_cs",
  "chem_157_12W_std_cs",
  "chem_140_12W_std_cs",
  "chem_38_12W_std_cs",
  "chem_127_12W_std_cs",
  "chem_3_12W_std_cs"
))

#Remove the week in the names
names(chems12)  <- sub('_12W_std_cs', '', names(chems12))

chems12$weekNum <- 12
chems12$weekID <- 1
chems12$patID <- 1:nrow(chems12)

#Rename total
names(chems12)[names(chems12) == "Total_12W"] <- "Total"

chems24 <- dplyr::select(OIF12and24,c(
  "chem_137_24W_std_cs",
  "chem_99_24W_std_cs",
  "chem_66_24W_std_cs",
  "chem_17_24W_std_cs",
  "chem_114_24W_std_cs",
  "chem_111_24W_std_cs",
  "chem_67_24W_std_cs",
  "chem_29_24W_std_cs",
  "chem_40_24W_std_cs",
  "chem_100_24W_std_cs",
  "chem_63_24W_std_cs",
  "chem_80_24W_std_cs",
  "chem_157_24W_std_cs",
  "chem_140_24W_std_cs",
  "chem_38_24W_std_cs",
  "chem_127_24W_std_cs",
  "chem_3_24W_std_cs"
))

names(chems24)  <- sub('_24W_std_cs', '', names(chems24))

chems24$weekNum <- 24
chems24$weekID <- 2
chems24$patID <- 1:nrow(chems24)

#Rename total
names(chems24)[names(chems24) == "Total_24W"] <- "Total"

#Merge
chemsLong <- rbind(chems12, chems24)

#Base Model

library(nlme)

mod<-lme(chem_3 ~ as.factor(weekID), random=list(patID=~1), data=chemsLong, na.action=na.omit)
sd1<-sqrt(getVarCov(mod)[1])
sd2<-mod$sigma
ICC<-sd1^2/(sd1^2+sd2^2)
ICC

#################################

#Bootstrapped 95% CIs

#It is necessary to enter each chemical individually

rm(list=setdiff(ls(), "chemsLong"))

set.seed(1)

numSamp <- 1000

capMat <- matrix(nrow = numSamp, ncol = 2)

for(p in 1:numSamp){
  
  nPat <- 20
  
  sample <- sample(c(1:nPat), replace = T)
  
  newDat <- matrix(ncol = ncol(chemsLong), nrow = nrow(chemsLong))
  
  for(i in 1:length(sample)){
    
    pair <- as.matrix(chemsLong[which(chemsLong$patID %in% sample[i]),])
    
    newDat[((2*i)-1):(2*i),] <- pair[1:2,]
    
  }
  
  #Conver to DF
  newDF <- as.data.frame(newDat)
  
  #Add nice names
  names(newDF) <- names(chemsLong)
  
  #Add new ID so each is treated as its own sample
  newDF$newID <- rep(1:nPat, each = 2)
  
  #The model and ICC
  mod<-lme(chem_3 ~ as.factor(weekID), random=list(newID=~1), data=newDF, na.action=na.omit)
  sd1<-sqrt(getVarCov(mod)[1])
  sd2<-mod$sigma
  ICC<-sd1^2/(sd1^2+sd2^2)
  
  capMat[p,1] <- p
  capMat[p,2] <- ICC
  
}

df <- as.data.frame(capMat)

hist(df$V2, breaks = 100)

mean(df$V2)
median(df$V2)
quantile(df$V2, probs = c(0.025, 0.975))
