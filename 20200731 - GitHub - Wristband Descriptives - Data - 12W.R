##############################

#Title: Assessment of Multipollutant Exposures during Pregnancy using Silicone Wristbands
#Authors: Brett T Doherty, John L Pearce, Kim A Anderson, Margaret R Karagas, Megan E Romano

##############################

#This file includes the code to prepare the data pertaining to the 12 gestational week wristbands

##############################

#Contents

#DDD1 - Create "chem" vars
#DDD2 - Create Chem Det Vars
#DDD3 - Create total detect count variable
#DDD4 - Create category-specific detect count variables

#DDD10 - Covariates

#DDD12 - Add batch-corrected chemicals
#DDD16 - Add batch-corrected chemicals, center & scaled

#####################

#Libraries

library(dplyr)

#####################

#Import Data

#Set directory
#XXX

#Import clean data
#XXX

#Import clean chemical class file
#XXX

#Set data to be cleaned, restrict to participants enrolled before 2019
toCleanData <- cleanData[which(cleanData$enroll_year %in% c(2017, 2018)),]

#Set classifications to be cleaned
toCleanClass <- cleanClass

##############################

#DDD1 - Create "chem" vars

#Create dataset that contains only the CASN variables
CASNVarsOld <- dplyr::select(toCleanData, starts_with("M_"))

#Sort dataset to ensure numeric/alphabetical order
CASNVarsOld<-CASNVarsOld[,order(colnames(CASNVarsOld),decreasing=FALSE)]

#Create column matrix of names of old CASN variables
CASNVarsOldNames <- as.matrix(names(CASNVarsOld))

#Create column matrix of names for new CASN variables
CASNVarsNewNames <- matrix(NA, nrow = length(CASNVarsOldNames), ncol = 1)
for(i in 1:length(CASNVarsOldNames)){CASNVarsNewNames[i] <- paste("chem", i, sep = "_")}

#Create matrix that includes both, to be used as an index
chemBook <- matrix(NA, nrow = length(CASNVarsOldNames), ncol = 2)
chemBook[,1] <- CASNVarsOldNames
chemBook[,2] <- CASNVarsNewNames

#Add these variables to the Classification
#First Sort
toCleanClass <- toCleanClass[order((toCleanClass$UNIQUE_ID__CASN),decreasing=FALSE),]
#Then add
toCleanClass$CASNOld <- chemBook[,1]
toCleanClass$chemName <- chemBook[,2]

#Apply new names to old dataset
CASNVarsNew <- CASNVarsOld
for(i in 1:ncol(CASNVarsNew)){colnames(CASNVarsNew)[i] <- CASNVarsNewNames[i]}

#Add these variables to the larger dataset
toCleanData <- cbind(toCleanData, CASNVarsNew)

#Write to Constructed - 2020
#XXX
#XXX

##############################

#DDD2 - Create Chem Det Vars

#Take new var dataset
CASNVarsNewDet <- CASNVarsNew

#Apply function to each column that creates a new variable: 1 if  chemX != 0, 0 if chemX = 0
CASNVarsNewDet[, paste0(names(CASNVarsNewDet),"_det")] <- lapply(CASNVarsNewDet, function(x) as.numeric(x!=0))

#Create dataset with these variables only
DetVars <- dplyr::select(CASNVarsNewDet, ends_with("_det"))

#Add these variables to the larger dataset
toCleanData <- cbind(toCleanData, DetVars)

#Write to Constructed - 2020
#XXX
#XXX

##############################

#DDD3 - Create total detect count variable

#Create new variable, Total, which is the sum of the det vars for each participant
DetVars$TotalDetect = rowSums(DetVars)

#Isolate this variable in preparation to append to larger dataset
TotalDet <- dplyr::select(DetVars, TotalDetect)

#Add variable to the larger dataset
toCleanData <- cbind(toCleanData, TotalDet)

#Write to Constructed - 2020
#XXX
#XXX

##############################

#DDD4 - Create category-specific detect count variables

#First, rename the variables in the class to make them easier to work with
names(toCleanClass)[names(toCleanClass) == 'Chemicals_in_Commerce'] <- 'CiC'
names(toCleanClass)[names(toCleanClass) == 'Consumer_Products'] <- 'CP'
names(toCleanClass)[names(toCleanClass) == 'Flame_Retardant'] <- 'FR'
names(toCleanClass)[names(toCleanClass) == 'Personal_Care'] <- 'PC'
names(toCleanClass)[names(toCleanClass) == 'Pesticides'] <- 'Pest'
names(toCleanClass)[names(toCleanClass) == 'Pharmacological'] <- 'Pharm'
names(toCleanClass)[names(toCleanClass) == 'Polycyclic_Aromatic_Hydrocarbon'] <- 'PAH'

#Make sets of each of the variables by their classifications
PAH <- toCleanClass[which(toCleanClass$PAH == 1),]
PAH_names <- as.matrix(PAH$chemName)

FR <- toCleanClass[which(toCleanClass$FR == 1),]
FR_names <- as.matrix(FR$chemName)

Pharm <- toCleanClass[which(toCleanClass$Pharm == 1),]
Pharm_names <- as.matrix(Pharm$chemName)

PC <- toCleanClass[which(toCleanClass$PC == 1),]
PC_names <- as.matrix(PC$chemName)

CiC <- toCleanClass[which(toCleanClass$CiC == 1),]
CiC_names <- as.matrix(CiC$chemName)

CP <- toCleanClass[which(toCleanClass$CP == 1),]
CP_names <- as.matrix(CP$chemName)

Pest <- toCleanClass[which(toCleanClass$Pest == 1),]
Pest_names <- as.matrix(Pest$chemName)

#Sum _det variables over each class;

#CiC
CiC_names_all <- matrix(NA, nrow = (length(CiC_names)), ncol = 1)
for(i in 1:length(CiC_names)){CiC_names_all[i,] <- paste(CiC_names[i],"_det", sep = "")}
CiCVarsDet <- toCleanData[,CiC_names_all]
CiCVarsDet$TotalCiC = rowSums(CiCVarsDet)
TotalCiC <- dplyr::select(CiCVarsDet, TotalCiC)

#CP
CP_names_all <- matrix(NA, nrow = (length(CP_names)), ncol = 1)
for(i in 1:length(CP_names)){CP_names_all[i,] <- paste(CP_names[i],"_det", sep = "")}
CPVarsDet  <- toCleanData[,CP_names_all]
CPVarsDet$TotalCP = rowSums(CPVarsDet)
TotalCP <- dplyr::select(CPVarsDet, TotalCP)

#FR
FR_names_all <- matrix(NA, nrow = (length(FR_names)), ncol = 1)
for(i in 1:length(FR_names)){FR_names_all[i,] <- paste(FR_names[i],"_det", sep = "")}
FRVarsDet  <- toCleanData[,FR_names_all]
FRVarsDet$TotalFR = rowSums(FRVarsDet)
TotalFR <- dplyr::select(FRVarsDet, TotalFR)

#PAH
PAH_names_all <- matrix(NA, nrow = (length(PAH_names)), ncol = 1)
for(i in 1:length(PAH_names)){PAH_names_all[i,] <- paste(PAH_names[i],"_det", sep = "")}
PAHVarsDet  <- toCleanData[,PAH_names_all]
PAHVarsDet$TotalPAH = rowSums(PAHVarsDet)
TotalPAH <- dplyr::select(PAHVarsDet, TotalPAH)

#PC
PC_names_all <- matrix(NA, nrow = (length(PC_names)), ncol = 1)
for(i in 1:length(PC_names)){PC_names_all[i,] <- paste(PC_names[i],"_det", sep = "")}
PCVarsDet  <- toCleanData[,PC_names_all]
PCVarsDet$TotalPC = rowSums(PCVarsDet)
TotalPC <- dplyr::select(PCVarsDet, TotalPC)

#Pest
Pest_names_all <- matrix(NA, nrow = (length(Pest_names)), ncol = 1)
for(i in 1:length(Pest_names)){Pest_names_all[i,] <- paste(Pest_names[i],"_det", sep = "")}
PestVarsDet <- toCleanData[,Pest_names_all]
PestVarsDet$TotalPest = rowSums(PestVarsDet)
TotalPest <- dplyr::select(PestVarsDet, TotalPest)

#Pharm
Pharm_names_all <- matrix(NA, nrow = (length(Pharm_names)), ncol = 1)
for(i in 1:length(Pharm_names)){Pharm_names_all[i,] <- paste(Pharm_names[i],"_det", sep = "")}
PharmVarsDet <- toCleanData[,Pharm_names_all]
PharmVarsDet$TotalPharm = rowSums(PharmVarsDet)
TotalPharm <- dplyr::select(PharmVarsDet, TotalPharm)

toCleanData <- cbind(toCleanData, TotalCiC, TotalCP, TotalFR, TotalPAH, TotalPC, TotalPest, TotalPharm)

#Write to Constructed - 2020
#XXX
#XXX

################################################

#DDD10 - Add new covariates

#Non-Hispanic White
toCleanData$mth_WNH <- ifelse(toCleanData$mth_race == "e. White" & toCleanData$mth_hisp == 0, 1, 0)
toCleanData$mth_WNHFlip <- ifelse(toCleanData$mth_race == "e. White" & toCleanData$mth_hisp == 0, 0, 1)

#2-level PNQRelatStat
toCleanData$PNQRelatStat2Lvl <- ifelse(toCleanData$PNQRelatStat == 1, 1, 0)
toCleanData$PNQRelatStat2LvlFlip <- ifelse(toCleanData$PNQRelatStat == 1, 0, 1)

#3-level PNQSchoolLvl
toCleanData$PNQSchoolLvl3Lvl <- ifelse(toCleanData$PNQSchoolLvl == 1 | toCleanData$PNQSchoolLvl == 2 | toCleanData$PNQSchoolLvl == 3, 0, toCleanData$PNQSchoolLvl)

#3-level PNQSchoolLvl
toCleanData$PNQSchoolLvl3LvlFlip <- ifelse(toCleanData$PNQSchoolLvl == 5, 0, 
                                           ifelse(toCleanData$PNQSchoolLvl == 4, 1, 2))

#Three 2-level PNQSchoolLVL
toCleanData$PNQSchoolLvl1 <- ifelse(toCleanData$PNQSchoolLvl == 1 | toCleanData$PNQSchoolLvl == 2 | toCleanData$PNQSchoolLvl == 3, 1, 
                              ifelse(toCleanData$PNQSchoolLvl == 4 | toCleanData$PNQSchoolLvl == 5, 0, NA))

toCleanData$PNQSchoolLvl2 <- ifelse(toCleanData$PNQSchoolLvl == 4, 1, 
                              ifelse(toCleanData$PNQSchoolLvl == 1 | toCleanData$PNQSchoolLvl == 2 | toCleanData$PNQSchoolLvl == 3 | toCleanData$PNQSchoolLvl == 5, 0, NA))

toCleanData$PNQSchoolLvl3 <- ifelse(toCleanData$PNQSchoolLvl == 5, 1, 
                              ifelse(toCleanData$PNQSchoolLvl == 1 | toCleanData$PNQSchoolLvl == 2 | toCleanData$PNQSchoolLvl == 3 | toCleanData$PNQSchoolLvl == 4, 0, NA))

#2-level parity
toCleanData$parity2Lvl <- ifelse(toCleanData$parity == 0, 1, 0)
toCleanData$parity2LvlFlip <- ifelse(toCleanData$parity == 0, 0, 1)

#4-level BMI
toCleanData$BMI4Lvl <- ifelse(toCleanData$maternalbmi == "NA", NA, 
                        ifelse(toCleanData$maternalbmi < 18.5, 0, 
                               ifelse(toCleanData$maternalbmi < 24.9, 1,  
                                      ifelse(toCleanData$maternalbmi < 29.9, 2, 3)))) 
#3-level BMI
toCleanData$BMI3Lvl <- ifelse(toCleanData$maternalbmi == "NA", NA, 
                                     ifelse(toCleanData$maternalbmi < 24.9, 0,  
                                            ifelse(toCleanData$maternalbmi < 29.9, 1, 2))) 

#WB Size 2-level
toCleanData$WBSize2Lvl <- ifelse(toCleanData$Wristband_Size == "regular", 1,
                              ifelse(toCleanData$Wristband_Size == "small", 0, NA))
toCleanData$WBSize2Lvl_flip <- ifelse(toCleanData$Wristband_Size == "regular", 0,
                                   ifelse(toCleanData$Wristband_Size == "small", 1, NA))

#Two-Level Smoking
toCleanData$smokegood <- ifelse(toCleanData$evercigpreg == 1, 1, 
                          ifelse(toCleanData$shsmokepreg == 1, 1, 
                                 ifelse(toCleanData$shsmokepreg == 0 & toCleanData$evercigpreg == 0, 0, NA)))
toCleanData$smokegoodFlip <- ifelse(toCleanData$evercigpreg == 1, 0, 
                                ifelse(toCleanData$shsmokepreg == 1, 0, 
                                       ifelse(toCleanData$shsmokepreg == 0 & toCleanData$evercigpreg == 0, 1, NA)))

#Gestage start and end in weeks
toCleanData$wb_start_gestage_weeks <- toCleanData$wb_start_gestage_days / 7
toCleanData$wb_end_gestage_weeks <- toCleanData$wb_end_gestage_days / 7

#Good Season
toCleanData$wb_start_season <- ifelse(toCleanData$wb_start_season == "", NA, 
                                         ifelse(toCleanData$wb_start_season == "fall", "fall", 
                                                ifelse(toCleanData$wb_start_season == "spring", "spring",
                                                       ifelse(toCleanData$wb_start_season == "summer", "summer",
                                                              ifelse(toCleanData$wb_start_season == "winter", "winter", NA)))))


#Season Indicators
toCleanData$seasonwearWint <- ifelse(is.na(toCleanData$wb_start_season), NA,
                                     ifelse(toCleanData$wb_start_season == "winter", 1, 0))
toCleanData$seasonwearSpri <- ifelse(is.na(toCleanData$wb_start_season), NA,
                                     ifelse(toCleanData$wb_start_season == "spring", 1, 0))
toCleanData$seasonwearSumm <- ifelse(is.na(toCleanData$wb_start_season), NA,
                                     ifelse(toCleanData$wb_start_season == "summer", 1, 0))
toCleanData$seasonwearFall <- ifelse(is.na(toCleanData$wb_start_season), NA,
                                     ifelse(toCleanData$wb_start_season == "fall", 1, 0))

#Difference between wristband and plasma sampling
toCleanData$diffPlasmaMinusWBStart <- toCleanData$blood_sample_gestage_days - toCleanData$wb_start_gestage_days

#Nail polish
toCleanData$nails <- ifelse(toCleanData$pnqnailfeett1 == 1 | toCleanData$pnqnailhandst1 == 1, 1,
                          ifelse(toCleanData$pnqnailfeett1 == 0 & toCleanData$pnqnailhandst1 == 0, 0, NA))

toCleanData$nailsFlip <- ifelse(toCleanData$pnqnailfeett1 == 1 | toCleanData$pnqnailhandst1 == 1, 0,
                            ifelse(toCleanData$pnqnailfeett1 == 0 & toCleanData$pnqnailhandst1 == 0, 1, NA))

#Handwashing
table(toCleanData$PNQWshHndsDay, useNA = "always")

#Replace impossible/implausible values with missing
toCleanData$PNQWshHndsDay <- ifelse(is.na(toCleanData$PNQWshHndsDay), NA, 
                                  ifelse(toCleanData$PNQWshHndsDay %in% c(-3,80,91515), NA, toCleanData$PNQWshHndsDay)) 

#Gardening and pesticides (full class, and specific pesiticides)
toCleanData$garden <- ifelse(toCleanData$PNQGardenT1 == 1, 1,
                           ifelse(toCleanData$PNQGardenT1 == 0, 0, NA))

toCleanData$gardenFlip <- ifelse(toCleanData$PNQGardenT1 == 1, 0,
                             ifelse(toCleanData$PNQGardenT1 == 0, 1, NA))

#Write to Constructed - 2020
#XXX
#XXX

################################################

#DDD12 - Add batch-corrected chemicals

#Plan:
  #Create a matrix of all the batch-chemical medians
  #Somehow cycle through the chemical matrix, deducting the appropriate value from each

#Step 1 - Create Matrix of batch-chem specific medians

chems <- dplyr::select(toCleanData, starts_with("chem"))
chems <- dplyr::select(chems, -contains("det"))
chems <- dplyr::select(chems, -contains("log"))

batch <- as.factor(toCleanData$Lab_Submission_Batch)

chemsBatch <- cbind(batch, chems)

capMat <- matrix(nrow = ncol(chems), ncol = 9)

for(i in 2:ncol(chemsBatch)){

capMat[i-1,] <- with(chemsBatch, tapply(chemsBatch[,i], chemsBatch[,1], median))

}

batchChemMeds <- capMat

#Step 2. Subtract the correct batch-chem median from each value.

table(toCleanData$Lab_Submission_Batch)

chems <- dplyr::select(toCleanData, starts_with("chem"))
chems <- dplyr::select(chems, -contains("det"))
chems <- dplyr::select(chems, -contains("log"))

batch <- as.factor(toCleanData$Lab_Submission_Batch)

chemsBatchOrig <- cbind(batch, chems)

chemsBatch <- cbind(chems,batch)

for(i in 1:(ncol(chemsBatch) - 1)){
  for(j in 1:nrow(chemsBatch)){
    
    chemsBatch[j,i] <- chemsBatch[j,i] - batchChemMeds[i, chemsBatch[j,203]]
    
  }
}

#remove the batch var
chemsBatch <- chemsBatch[,1:202]

#Rename the variables
colnames(chemsBatch) <- paste(colnames(chemsBatch), "_std", sep = "")

#Add to toCleanData
toCleanData <- cbind(toCleanData, chemsBatch)

#Write to Constructed - 2020
#XXX
#XXX

###################################################################

#DDD16 - Add batch-corrected, standardized variables

batchStan <- dplyr::select(toCleanData, ends_with("std"))
batchStan <- dplyr::select(batchStan, -contains("log2"))

batchStanStan <- as.data.frame(scale(batchStan, center = T, scale = T))

#Rename the variables
colnames(batchStanStan) <- paste(colnames(batchStanStan), "_cs", sep = "")

#Add to toCleanData
toCleanData <- cbind(toCleanData, batchStanStan)

#Write to Constructed - 2020
#XXX
#XXX

