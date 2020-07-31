##############################

#Title: Assessment of Multipollutant Exposures during Pregnancy using Silicone Wristbands
#Authors: Brett T Doherty, John L Pearce, Kim A Anderson, Margaret R Karagas, Megan E Romano

##############################

#This file includes the code to analyze the data pertaining to the 12 gestational week wristbands

##############################

#Contents

#AAA - Population

#BBB - Chemical Concentrations

#CCC - Chemical Category Summary Statistics

#DDD - Correlations

#EEE - Multivariable Linear Regression

#FFF - Specific Exposure-Outcome Associations

#GGG - SOM

#HHH - SOM and Chemical Concentrations and Covariates

#III - Modeling the pollutant concentrations as a function of season

######################################

#Libraries

library(dplyr)

######################################

#Data

#Set directory
#XXX

#Read cleaned data
#XXX

#Read cleaned class data
#XXX

######################################

#AAA - Population

library(janitor)

summary(constData$enrollment_age)
summary(constData$maternalbmi)
tabyl(constData$PNQSchoolLvl3Lvl)
tabyl(constData$PNQRelatStat2LvlFlip)
tabyl(constData$mth_WNHFlip)
tabyl(constData$parity2LvlFlip)
tabyl(constData$smokegood)

summary(constData$wb_start_gestage_weeks)
summary(constData$Days_worn)
tabyl(constData$Wristband_Size)
tabyl(constData$wb_start_season)

##############################

#BBB - Chemical Concentrations

#Summary detection (i.e., % detected for each chemical)

chems <- select(constData, contains("chem"))
chemsdet <- select(chems, contains("det"))
chemsdet <- chemsdet[complete.cases(chemsdet),]

DET<-function (x){length(which(x == 1))}
DETs<-round(apply(chemsdet,2,FUN=DET),2)
DETs.per <- round(100*(DETs/dim(chemsdet)[1]),2)

#Basic descriptives

chems <- select(constData, contains("chem"))
chemsconc <- select(chems, -contains("det"))
chemsconc <- select(chemsconc, -contains("log"))
chemsconc <- select(chemsconc, -contains("std"))
chemsconc <- chemsconc[complete.cases(chemsconc),]

p50 <- round(apply(chemsconc,2,FUN='quantile',probs = 0.50, na.rm=TRUE),0)
p25 <- round(apply(chemsconc,2,FUN='quantile',probs = 0.25, na.rm=TRUE),0)
p75 <- round(apply(chemsconc,2,FUN='quantile',probs = 0.75, na.rm=TRUE),0)
p95 <- round(apply(chemsconc,2,FUN='quantile',probs = 0.95, na.rm=TRUE),0)
max <- round(apply(chemsconc,2,FUN='max', na.rm=TRUE),0)

chemsstats <- cbind(p50, p25, p75, p95, max)

##############################

#CCC - Chemical Category Summary Statistics

capMat <- matrix(nrow = 8, ncol = 6)

a <- summary(constData$TotalDetect)
capMat[1,1] <- "Total"
capMat[1,2] <- a[3]
capMat[1,3] <- a[2]
capMat[1,4] <- a[5]
capMat[1,5] <- a[1]
capMat[1,6] <- a[6]

a <- summary(constData$TotalCiC)
capMat[2,1] <- "CiC"
capMat[2,2] <- a[3]
capMat[2,3] <- a[2]
capMat[2,4] <- a[5]
capMat[2,5] <- a[1]
capMat[2,6] <- a[6]

a <- summary(constData$TotalPC)
capMat[3,1] <- "PC"
capMat[3,2] <- a[3]
capMat[3,3] <- a[2]
capMat[3,4] <- a[5]
capMat[3,5] <- a[1]
capMat[3,6] <- a[6]

a <- summary(constData$TotalPest)
capMat[4,1] <- "Pest"
capMat[4,2] <- a[3]
capMat[4,3] <- a[2]
capMat[4,4] <- a[5]
capMat[4,5] <- a[1]
capMat[4,6] <- a[6]

a <- summary(constData$TotalFR)
capMat[5,1] <- "FR"
capMat[5,2] <- a[3]
capMat[5,3] <- a[2]
capMat[5,4] <- a[5]
capMat[5,5] <- a[1]
capMat[5,6] <- a[6]

a <- summary(constData$TotalPAH)
capMat[6,1] <- "PAH"
capMat[6,2] <- a[3]
capMat[6,3] <- a[2]
capMat[6,4] <- a[5]
capMat[6,5] <- a[1]
capMat[6,6] <- a[6]

a <- summary(constData$TotalCP)
capMat[7,1] <- "CP"
capMat[7,2] <- a[3]
capMat[7,3] <- a[2]
capMat[7,4] <- a[5]
capMat[7,5] <- a[1]
capMat[7,6] <- a[6]

a <- summary(constData$TotalPharm)
capMat[8,1] <- "Pharm"
capMat[8,2] <- a[3]
capMat[8,3] <- a[2]
capMat[8,4] <- a[5]
capMat[8,5] <- a[1]
capMat[8,6] <- a[6]

##############################

#DDD - Correlations among >60%

chems <- select(constData, c(
  "chem_178_std_cs",
  "chem_127_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_146_std_cs",
  "chem_142_std_cs",
  "chem_86_std_cs",
  "chem_36_std_cs",
  "chem_52_std_cs",
  "chem_80_std_cs",
  "chem_128_std_cs",
  "chem_102_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs",
  "chem_50_std_cs",
  "chem_162_std_cs"
  ))

a <- cor(chems, method = "spearman", use = "pairwise.complete.obs")
a[lower.tri(a)] <- ""

##############################

#Heatmap Figure

chems60 <- select(constData, c(
  "chem_178_std_cs",
  "chem_127_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_146_std_cs",
  "chem_142_std_cs",
  "chem_86_std_cs",
  "chem_36_std_cs",
  "chem_52_std_cs",
  "chem_80_std_cs",
  "chem_128_std_cs",
  "chem_102_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs",
  "chem_50_std_cs",
  "chem_162_std_cs"
))

colnames(chems60) <- c(
  "Di-n-butyl phthalate",
  "Galaxolide",
  "Diisobutyl phthalate",
  "Butyl benzyl phthalate",
  "Lilial",
  "Benzyl salicylate",
  "Tonalide",
  "N,N-Diethyl-m-toluamide",
  "Benzophenone",
  "Benzyl benzoate",
  "Ethylene brassylate",
  "Di-n-nonyl phthalate",
  "Permethrin",
  "Diethyl phthalate",
  "Butylated hydroxyanisole",
  "2,4-Di-tert-butylphenol"
)

chem60cor <- cor(chems60, method = "spearman")

library(gplots)

chem60corR <- round(chem60cor, 2)

heatmap.2(chem60corR, 
          cellnote = chem60corR,
          col = bluered(100),
          notecol = "black",
          trace = "none",
          density.info = "none",
          dendrogram = "none",
          margins = c(16,16),
          key = TRUE,
          keysize = 1,
          key.title = "",
          key.xlab = "Spearman Correlation",
          cexRow = 1.5,
          cexCol = 1.5)

##############################

#EEE - Multivariable Linear Regression

mvlr <- dplyr::select(constData, c(
  "TotalDetect",
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3Lvl",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2LvlFlip",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Convert to factors and numeric
mvlr$TotalDetect <- as.numeric(mvlr$TotalDetect)
mvlr$enrollment_age <- as.numeric(mvlr$enrollment_age)
mvlr$maternalbmi <- as.numeric(mvlr$maternalbmi)
mvlr$PNQSchoolLvl3Lvl <- as.factor(mvlr$PNQSchoolLvl3Lvl)
mvlr$PNQRelatStat2LvlFlip <- as.factor(mvlr$PNQRelatStat2LvlFlip)
mvlr$mth_WNHFlip <- as.factor(mvlr$mth_WNHFlip)
mvlr$parity2LvlFlip <- as.factor(mvlr$parity2LvlFlip)
mvlr$smokegood <- as.factor(mvlr$smokegood)
mvlr$wb_start_gestage_weeks <- as.numeric(mvlr$wb_start_gestage_weeks)
mvlr$wb_start_season <- as.factor(mvlr$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(mvlr, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove TotalDetect from the variables used to impute the others
predM[, c("TotalDetect")] <- 0

#Run the Imputation
set.seed(1)
mvlrVarsImp <- mice(mvlr, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
fitMI <- with(mvlrVarsImp, lm(TotalDetect ~ enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3Lvl) +
                                 as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2LvlFlip) +
                                 as.factor(smokegood) + wb_start_gestage_weeks + as.factor(wb_start_season)))

pooled <- pool(fitMI)

out <- summary(pooled)

##############################

#FFF - Specific Exposure-Outcome Associations

##############################

#Nail polish in first 3 months - TPP + phthalates

table(constData$nails, useNA = "always")

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "nails",
  
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
#TotalPC
# TPP	                    chem_3_std_cs
# Di-n-butyl phthalate	  chem_178_std_cs
# Diisobutyl phthalate	  chem_85_std_cs
# Butyl benzyl phthalate	chem_20_std_cs
# Di-n-nonyl phthalate	  chem_102_std_cs
# Diethyl phthalate	      chem_181_std_cs

chems <- dplyr::select(constData, c(
  "TotalPC",
  "chem_3_std_cs",
  "chem_178_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_102_std_cs",
  "chem_181_std_cs"
  ))

#Combine
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$TotalPC <- as.numeric(constDataVars$TotalPC)
constDataVars$chem_3_std_cs <- as.numeric(constDataVars$chem_3_std_cs)
constDataVars$chem_178_std_cs <- as.numeric(constDataVars$chem_178_std_cs)
constDataVars$chem_85_std_cs <- as.numeric(constDataVars$chem_85_std_cs)
constDataVars$chem_20_std_cs <- as.numeric(constDataVars$chem_20_std_cs)
constDataVars$chem_102_std_cs <- as.numeric(constDataVars$chem_102_std_cs)
constDataVars$chem_181_std_cs <- as.numeric(constDataVars$chem_181_std_cs)

constDataVars$nails <- as.factor(constDataVars$nails)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Chemical Outcome from the variables used to impute the others
predM[, c("TotalPC")] <- 0
predM[, c("chem_3_std_cs")] <- 0
predM[, c("chem_178_std_cs")] <- 0
predM[, c("chem_85_std_cs")] <- 0
predM[, c("chem_20_std_cs")] <- 0
predM[, c("chem_102_std_cs")] <- 0
predM[, c("chem_181_std_cs")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
#TotalPC
# TPP	                    chem_3_std_cs
# Di-n-butyl phthalate	  chem_178_std_cs
# Diisobutyl phthalate	  chem_85_std_cs
# Butyl benzyl phthalate	chem_20_std_cs
# Di-n-nonyl phthalate	  chem_102_std_cs
# Diethyl phthalate	      chem_181_std_cs

#TotalPC
fitMI <- with(constDataVarsImp, lm(TotalPC ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# TPP	                    chem_3_std_cs
fitMI <- with(constDataVarsImp, lm(chem_3_std_cs ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Di-n-butyl phthalate	  chem_178_std_cs
fitMI <- with(constDataVarsImp, lm(chem_178_std_cs ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                 as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                   wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Diisobutyl phthalate	  chem_85_std_cs
fitMI <- with(constDataVarsImp, lm(chem_85_std_cs ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Butyl benzyl phthalate	chem_20_std_cs
fitMI <- with(constDataVarsImp, lm(chem_20_std_cs ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Di-n-nonyl phthalate	  chem_102_std_cs
fitMI <- with(constDataVarsImp, lm(chem_102_std_cs ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Diethyl phthalate	      chem_181_std_cs
fitMI <- with(constDataVarsImp, lm(chem_181_std_cs ~ as.factor(nails) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

##############################

#Handwashing and total # chemicals

table(constData$PNQWshHndsDay, useNA = "always")
summary(constData$PNQWshHndsDay)

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "TotalDetect",
  
  "PNQWshHndsDay",
  
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Convert to factors and numeric
constDataVars$TotalDetect <- as.numeric(constDataVars$TotalDetect)

constDataVars$PNQWshHndsDay <- as.numeric(constDataVars$PNQWshHndsDay)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove variables used to impute the others
predM[, c("TotalDetect")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
fitMI <- with(constDataVarsImp, lm(TotalDetect ~ PNQWshHndsDay + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

##############################

#Gardening and pesticides (full class, and specific pesiticides)

table(constData$garden, useNA = "always")

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "garden",
  
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
#TotalPest
# Di-n-butyl phthalate	      chem_178_std_cs
# N,N-Diethyl-m-toluamide	    chem_36_std_cs
# Benzyl benzoate	            chem_80_std_cs
# Permethrin	                chem_201_std_cs
# Diethyl phthalate	          chem_181_std_cs

chems <- dplyr::select(constData, c(
  "TotalPest",
  "chem_178_std_cs",
  "chem_36_std_cs",
  "chem_80_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs"
))

#Add back
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$TotalPest <- as.numeric(constDataVars$TotalPest)
constDataVars$chem_178_std_cs <- as.numeric(constDataVars$chem_178_std_cs)
constDataVars$chem_36_std_cs <- as.numeric(constDataVars$chem_36_std_cs)
constDataVars$chem_80_std_cs <- as.numeric(constDataVars$chem_80_std_cs)
constDataVars$chem_201_std_cs <- as.numeric(constDataVars$chem_201_std_cs)
constDataVars$chem_181_std_cs <- as.numeric(constDataVars$chem_181_std_cs)

constDataVars$garden <- as.factor(constDataVars$garden)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Outcome from the variables used to impute the others
#TotalPest
# Di-n-butyl phthalate	      chem_178_std_cs
# N,N-Diethyl-m-toluamide	    chem_36_std_cs
# Benzyl benzoate	            chem_80_std_cs
# Permethrin	                chem_201_std_cs
# Diethyl phthalate	          chem_181_std_cs
predM[, c("TotalPest")] <- 0
predM[, c("chem_178_std_cs")] <- 0
predM[, c("chem_36_std_cs")] <- 0
predM[, c("chem_80_std_cs")] <- 0
predM[, c("chem_201_std_cs")] <- 0
predM[, c("chem_181_std_cs")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
#TotalPest
# Di-n-butyl phthalate	      chem_178_std_cs
# N,N-Diethyl-m-toluamide	    chem_36_std_cs
# Benzyl benzoate	            chem_80_std_cs
# Permethrin	                chem_201_std_cs
# Diethyl phthalate	          chem_181_std_cs

#TotalPest
fitMI <- with(constDataVarsImp, lm(TotalPest ~ as.factor(garden) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                 as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                 wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Di-n-butyl phthalate	      chem_178_std_cs
fitMI <- with(constDataVarsImp, lm(chem_178_std_cs ~ as.factor(garden) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# N,N-Diethyl-m-toluamide	    chem_36_std_cs
fitMI <- with(constDataVarsImp, lm(chem_36_std_cs ~ as.factor(garden) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Benzyl benzoate	            chem_80_std_cs
fitMI <- with(constDataVarsImp, lm(chem_80_std_cs ~ as.factor(garden) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Permethrin	                chem_201_std_cs
fitMI <- with(constDataVarsImp, lm(chem_201_std_cs ~ as.factor(garden) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

# Diethyl phthalate	          chem_181_std_cs
fitMI <- with(constDataVarsImp, lm(chem_181_std_cs ~ as.factor(garden) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season)))
pooled <- pool(fitMI)
summary(pooled)

##############################

#Season + benzophenone

table(constData$wb_start_season, useNA = "always")

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
chems <- dplyr::select(constData, c(
  "chem_52_std_cs"))

#Add back
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$chem_52_std_cs <- as.numeric(constDataVars$chem_52_std_cs)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Outcome from the variables used to impute the others
predM[, c("chem_52_std_cs")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
fitMI <- with(constDataVarsImp, lm(chem_52_std_cs ~ as.factor(wb_start_season) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks))
pooled <- pool(fitMI)
summary(pooled)

##############################

#Season + Total Pesticide

table(constData$wb_start_season, useNA = "always")

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
chems <- dplyr::select(constData, c(
  "TotalPest"))

#Add back
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$TotalPest <- as.numeric(constDataVars$TotalPest)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Outcome from the variables used to impute the others
predM[, c("TotalPest")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
fitMI <- with(constDataVarsImp, lm(TotalPest ~ as.factor(wb_start_season) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks))
pooled <- pool(fitMI)
summary(pooled)

##############################

#Season + N,N-Diethyl-m-toluamide

table(constData$wb_start_season, useNA = "always")

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
chems <- dplyr::select(constData, c(
  "chem_36_std_cs"))

#Add back
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$chem_36_std_cs <- as.numeric(constDataVars$chem_36_std_cs)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Outcome from the variables used to impute the others
predM[, c("chem_36_std_cs")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
fitMI <- with(constDataVarsImp, lm(chem_36_std_cs ~ as.factor(wb_start_season) + enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2Lvl) + as.factor(smokegood) +
                                     wb_start_gestage_weeks))
pooled <- pool(fitMI)
summary(pooled)


##############################

#Parity + TPHP

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2LvlFlip",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
chems <- dplyr::select(constData, c(
  "chem_3_std_cs"))

#Add back
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$chem_3_std_cs <- as.numeric(constDataVars$chem_3_std_cs)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2LvlFlip <- as.factor(constDataVars$parity2LvlFlip)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Outcome from the variables used to impute the others
predM[, c("chem_3_std_cs")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The Multiple Imputed Linear Model
fitMI <- with(constDataVarsImp, lm(chem_3_std_cs ~ enrollment_age + maternalbmi + as.factor(PNQSchoolLvl3LvlFlip) +
                                     as.factor(PNQRelatStat2LvlFlip) + as.factor(mth_WNHFlip) + as.factor(parity2LvlFlip) + as.factor(smokegood) +
                                     wb_start_gestage_weeks + as.factor(wb_start_season) ))
pooled <- pool(fitMI)
summary(pooled)

#######################################################

#GGG - SOM - Chemicals

#XXX

#XXX

#XXX

#constClassSOM
  #Dataset where I add a variable "oneClass" to constClass, where the "primary" class of the chemicals is indicated
  #I also add 'Detect' to this
#XXX

#SOM Code Functions
#XXX

chems <- dplyr::select(constData, c(
  "chem_178_std_cs",
  "chem_127_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_146_std_cs",
  "chem_142_std_cs",
  "chem_86_std_cs",
  "chem_36_std_cs",
  "chem_52_std_cs",
  "chem_80_std_cs",
  "chem_128_std_cs",
  "chem_102_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs",
  "chem_50_std_cs",
  "chem_162_std_cs"
  ))

#Standardize
chems <-scale(chems, center=TRUE, scale =TRUE)

#Back to Dataframe
chems <- as.data.frame(chems)

#Look at data with histograms
chemNames <- c(
  "Di-n-butyl phthalate",
  "Galaxolide",
  "Diisobutyl phthalate",
  "Butyl benzyl phthalate",
  "Lilial",
  "Benzyl salicylate",
  "Tonalide",
  "N,N-Diethyl-m-toluamide",
  "Benzophenone",
  "Benzyl benzoate",
  "Ethylene brassylate",
  "Di-n-nonyl phthalate",
  "Permethrin",
  "Diethyl phthalate",
  "Butylated hydroxyanisole",
  "2,4-Di-tert-butylphenol"
)

#Nice Name
SOMdata <- chems

#Add names
names(SOMdata) <- chemNames

#Prep class in preparation for the figure
class60 <- constClassSOM[which(constClassSOM$Detect > 60),]
class60 <- class60[order(class60$oneClass, class60$Compound_Name),]

#Get vector of the chemical names straight from chemInfo25
a <- as.data.frame(class60$Compound_Name)
#XXX

#Reorder accordingly
SOMdata <- dplyr::select(SOMdata, c(
  "Butyl benzyl phthalate",
  "Di-n-butyl phthalate",
  "Di-n-nonyl phthalate",
  "Diethyl phthalate",
  "Diisobutyl phthalate",
  "Benzyl benzoate",
  "N,N-Diethyl-m-toluamide",
  "Permethrin",
  "2,4-Di-tert-butylphenol",
  "Galaxolide",
  "Benzophenone",
  "Benzyl salicylate",
  "Butylated hydroxyanisole",
  "Ethylene brassylate",
  "Lilial",
  "Tonalide"
))

################################

#SOM evaluation functions

#Var eval
var.eval(SOMdata)

#Check correlation structure
cor.str.eval(SOMdata)

#Check grouping structure
grp.str.eval(SOMdata, kmn=2, kmx=20)

#map evaluations
map.size.eval(SOMdata, 7)

##############################

#Make SOM

data <- as.matrix(SOMdata)

somx=3
somy=4
ksize=somx*somy

set.seed(2)
fin.som<-som(data, grid=somgrid(xdim=somx,ydim=somy, "rectangular"), 
             rlen=50000, alpha=c(0.05,0.01)) 

som.summ <- som.fit.summ(fin.som)
som.summ$SOM_N

##############################

#Star Figure

par(mfrow=c(1,1), mar=c(5,4,2,2), family="serif")
profiles<-fin.som$codes[[1]]
IDs=1:dim(profiles)[1]
p=dim(profiles)[2]
var.labs<-dimnames(profiles)[[2]]

#Number of colors needed for each class colors
#5, 3, 2, 6

library(RColorBrewer)
col5 <- brewer.pal(n = 5, name = "Blues")
col3 <- brewer.pal(n = 3, name = "Greens")
col2 <- c("#CCCCCC", "#333333")
col6 <- brewer.pal(n = 6, name = "Reds")

cols <- c(col5, col3, col2, col6)

palette(cols)

stars(profiles, locations=fin.som$grid$pts, draw.segments=TRUE, axes=FALSE, scale=TRUE, 
      len = 0.35, labels=NULL,
      ylim=c(0.25,somy+0.5), xlim=c(0.75,somx+1.75))

symbols(fin.som$grid$pts[, 1], fin.som$grid$pts[, 2],
        circles = rep(.48, nrow(fin.som$grid$pts)),
        inches = FALSE, add = TRUE,
        fg = "black", bg = NA)

#Reorder IDs per reviewer request
newIDs <- c(10, 11, 12, 7, 8, 9, 4, 5, 6, 1, 2, 3)

text(x=fin.som$grid$pts[, 1], y=fin.som$grid$pts[, 2]+.40,
     labels=paste("[", newIDs ,"]", sep=""), font=2, cex=1)

legend("right", legend=var.labs, pch=15, col = palette(),
       pt.cex=1.5,cex=1.35, angle=45, ncol=1, inset=0.08, horiz=F, y.intersp = 0.40, bty = "n")

text(x=fin.som$grid$pts[, 1], y=fin.som$grid$pts[, 2] - 0.40,
     labels=paste("(n = ", som.summ$SOM_N ,")", sep=""), font=2, cex=1)

##############################

#HHH - SOM and Chemical Concentrations and Covariates

#The ID Variable
som.summ$SOM_CLASSIF$ID

#Add to dataset
SOMdata$SOMID <- som.summ$SOM_CLASSIF$ID

#Confirm
table(SOMdata$SOMID)

#Chemicals
capMat <- matrix(ncol = somx*somy, nrow = (ncol(SOMdata) - 1))

for(i in 1:(ncol(SOMdata) - 1)){
  capMat[i,] <- tapply(SOMdata[,i], SOMdata$SOMID, median, na.rm = T)
}

#XXX

###################################

#Covariates

#Add to dataset
constData$SOMID <- som.summ$SOM_CLASSIF$ID

attach(constData)

b <- tapply(enrollment_age, SOMID, median, na.rm = T)
c <- tapply(maternalbmi, SOMID, median, na.rm = T)
d <- tapply(PNQSchoolLvl1, SOMID, mean, na.rm = T)
e <- tapply(PNQSchoolLvl2, SOMID, mean, na.rm = T)
f <- tapply(PNQSchoolLvl3, SOMID, mean, na.rm = T)
g <- tapply(PNQRelatStat2Lvl, SOMID, mean, na.rm = T)
h <- tapply(PNQRelatStat2LvlFlip, SOMID, mean, na.rm = T)
i <- tapply(mth_WNH, SOMID, mean, na.rm = T)
j <- tapply(mth_WNHFlip, SOMID, mean, na.rm = T)
k <- tapply(parity2Lvl, SOMID, mean, na.rm = T)
l <- tapply(parity2LvlFlip, SOMID, mean, na.rm = T)
m <- tapply(smokegoodFlip, SOMID, mean, na.rm = T)
n <- tapply(smokegood, SOMID, mean, na.rm = T)
o <- tapply(wb_start_gestage_weeks, SOMID, median, na.rm = T)
s <- tapply(seasonwearWint, SOMID, mean, na.rm = T)
t <- tapply(seasonwearSpri, SOMID, mean, na.rm = T)
u <- tapply(seasonwearSumm, SOMID, mean, na.rm = T)
v <- tapply(seasonwearFall, SOMID, mean, na.rm = T)
w <- tapply(nails, SOMID, mean, na.rm = T)
x <- tapply(nailsFlip, SOMID, mean, na.rm = T)
y <- tapply(PNQWshHndsDay, SOMID, median, na.rm = T)
z <- tapply(garden, SOMID, mean, na.rm = T)
aa <- tapply(gardenFlip, SOMID, mean, na.rm = T)

SOMChemMeds <- rbind(b,c,d,e,f,g,h,i,j,k,l,m,n,o,s,t,u,v,w,x,y,z,aa)

#XXX

##############################

#III - Modeling the pollutant concentrations as a function of season

##############################

#Select covariates
constDataVars <- dplyr::select(constData, c(
  "enrollment_age",
  "maternalbmi",
  "PNQSchoolLvl3LvlFlip",
  "PNQRelatStat2LvlFlip",
  "mth_WNHFlip",
  "parity2Lvl",
  "smokegood",
  "wb_start_gestage_weeks",
  "wb_start_season"
))

#Select chemicals
chems <- dplyr::select(constData, c(
  "chem_178_std_cs",
  "chem_127_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_146_std_cs",
  "chem_142_std_cs",
  "chem_86_std_cs",
  "chem_36_std_cs",
  "chem_52_std_cs",
  "chem_80_std_cs",
  "chem_128_std_cs",
  "chem_102_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs",
  "chem_50_std_cs",
  "chem_162_std_cs"))

#Add back
constDataVars <- cbind(chems, constDataVars)

#Convert to factors and numeric
constDataVars$chem_178_std_cs <- as.numeric(constDataVars$chem_178_std_cs)
constDataVars$chem_127_std_cs <- as.numeric(constDataVars$chem_127_std_cs)
constDataVars$chem_85_std_cs <- as.numeric(constDataVars$chem_85_std_cs)
constDataVars$chem_20_std_cs <- as.numeric(constDataVars$chem_20_std_cs)
constDataVars$chem_146_std_cs <- as.numeric(constDataVars$chem_146_std_cs)
constDataVars$chem_142_std_cs <- as.numeric(constDataVars$chem_142_std_cs)
constDataVars$chem_86_std_cs <- as.numeric(constDataVars$chem_86_std_cs)
constDataVars$chem_36_std_cs <- as.numeric(constDataVars$chem_36_std_cs)
constDataVars$chem_52_std_cs <- as.numeric(constDataVars$chem_52_std_cs)
constDataVars$chem_80_std_cs <- as.numeric(constDataVars$chem_80_std_cs)
constDataVars$chem_128_std_cs <- as.numeric(constDataVars$chem_128_std_cs)
constDataVars$chem_102_std_cs <- as.numeric(constDataVars$chem_102_std_cs)
constDataVars$chem_201_std_cs <- as.numeric(constDataVars$chem_201_std_cs)
constDataVars$chem_181_std_cs <- as.numeric(constDataVars$chem_181_std_cs)
constDataVars$chem_50_std_cs <- as.numeric(constDataVars$chem_50_std_cs)
constDataVars$chem_162_std_cs <- as.numeric(constDataVars$chem_162_std_cs)

constDataVars$enrollment_age <- as.numeric(constDataVars$enrollment_age)
constDataVars$maternalbmi <- as.numeric(constDataVars$maternalbmi)
constDataVars$PNQSchoolLvl3LvlFlip <- as.factor(constDataVars$PNQSchoolLvl3LvlFlip)
constDataVars$PNQRelatStat2LvlFlip <- as.factor(constDataVars$PNQRelatStat2LvlFlip)
constDataVars$mth_WNHFlip <- as.factor(constDataVars$mth_WNHFlip)
constDataVars$parity2Lvl <- as.factor(constDataVars$parity2Lvl)
constDataVars$smokegood <- as.factor(constDataVars$smokegood)
constDataVars$wb_start_gestage_weeks <- as.numeric(constDataVars$wb_start_gestage_weeks)
constDataVars$wb_start_season <- as.factor(constDataVars$wb_start_season)

#Impute with mice
library(mice)

#Initialize
init <- mice(constDataVars, maxit = 0)
meth <- init$method
predM <- init$predictorMatrix

#Remove Outcome from the variables used to impute the others
predM[, c(
  "chem_178_std_cs",
  "chem_127_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_146_std_cs",
  "chem_142_std_cs",
  "chem_86_std_cs",
  "chem_36_std_cs",
  "chem_52_std_cs",
  "chem_80_std_cs",
  "chem_128_std_cs",
  "chem_102_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs",
  "chem_50_std_cs",
  "chem_162_std_cs")] <- 0

#Run the Imputation
set.seed(1)
constDataVarsImp <- mice(constDataVars, method = meth, predictorMatrix = predM, m = 25)

#The loop

chems <- c(
  "chem_178_std_cs",
  "chem_127_std_cs",
  "chem_85_std_cs",
  "chem_20_std_cs",
  "chem_146_std_cs",
  "chem_142_std_cs",
  "chem_86_std_cs",
  "chem_36_std_cs",
  "chem_52_std_cs",
  "chem_80_std_cs",
  "chem_128_std_cs",
  "chem_102_std_cs",
  "chem_201_std_cs",
  "chem_181_std_cs",
  "chem_50_std_cs",
  "chem_162_std_cs")

capMat <- matrix(nrow = length(chems), ncol = 10)

colnames(capMat) <- c("Chemical", 
                      "Spring_Beta", "Spring_95L", "Spring_95U",
                      "Summer_Beta", "Summer_95L", "Summer_95U",
                      "Winter_Beta", "Winter_95L", "Winter_95U")

covars <- c("as.factor(wb_start_season)", "enrollment_age", "maternalbmi", "as.factor(PNQSchoolLvl3LvlFlip)",
            "as.factor(PNQRelatStat2LvlFlip)", "as.factor(mth_WNHFlip)", "as.factor(parity2Lvl)", "as.factor(smokegood)",
            "wb_start_gestage_weeks")

for(i in 1:length(chems)){
  
  fitMI <- with(constDataVarsImp, lm(as.formula(paste(chems[i], paste(covars, collapse=" + "), sep=" ~ "))))
  
  pooled <- pool(fitMI)
  b <- summary(pooled)
  
  capMat[i,1] <- chems[i]
  
  capMat[i,2] <- b$estimate[2]
  capMat[i,3] <- b$estimate[2] - b$std.error[2]*1.96
  capMat[i,4] <- b$estimate[2] + b$std.error[2]*1.96
  
  capMat[i,5] <- b$estimate[3]
  capMat[i,6] <- b$estimate[3] - b$std.error[3]*1.96
  capMat[i,7] <- b$estimate[3] + b$std.error[3]*1.96
  
  capMat[i,8] <- b$estimate[4]
  capMat[i,9] <- b$estimate[4] - b$std.error[4]*1.96
  capMat[i,10] <- b$estimate[4] + b$std.error[4]*1.96
  
}

#XXX
