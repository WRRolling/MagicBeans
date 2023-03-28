# Results:Genomic prediction using US PI panels to generate or validate models
# Author: William Rolling
# Corresponding Author: Leah McHale (mchale.21@osu.edu)
#####################################################################################
# Section 1.0: Prepare R environment
# 1.1 Clear R environment
rm(list=ls())
# 1.2 List packages need to complete analysis
Package.List <- c("BGLR",
                  "bigmemory",
                  "data.table",
                  "bigmemory",
                  "biganalytics",
                  "tidyverse",
                  "parallel",
                  "MASS")
# MASS and dplyr select do not play well together
# 1.3 Install and/or load R packages
# Script from: https://vbaliga.github.io
package.check <- lapply(
  Package.List ,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})
# 1.4 Setworking directory
# example for someone who downloads data to folder in desktop <---- Set pathway to folder$    setwd("~/../Desktop/Supplemental Data 1/")
#####################################################################################
# Section 2.0: Load genotypic data into R environment
# 2.1 Load genotype files for three experiments
SK.2016 <- read_csv("../Input.data/SK.2016.hmp.csv", col_names=TRUE)
C2.2018 <- read_csv("../Input.data/C2.2018.hmp.csv", col_names=TRUE)
OH.2018 <- read_csv("../Input.data/OH.2018.hmp.csv", col_names=TRUE)
# 2.2 rename first column to match SK.2016 column name
colnames(C2.2018)[1] <- "rs"
colnames(OH.2018)[1] <- "rs"
# 2.3 Match SNPs across all three populations in SK.2016 Experiment
SK.2016.formatted <- filter(SK.2016, rs %in% C2.2018$rs)
SK.2016.formatted <- filter(SK.2016.formatted, rs %in% OH.2018$rs)
SK.2016.formatted <- SK.2016.formatted[order(SK.2016.formatted$rs),]
# 2.4 Match SNPs across all three populations in OH.2018 and C2.2018
C2.2018.formatted <- filter(C2.2018, rs %in% SK.2016.formatted$rs)
C2.2018.formatted <- C2.2018.formatted[order(C2.2018.formatted$rs),]
OH.2018.formatted <- filter(OH.2018,rs %in% SK.2016.formatted$rs)
OH.2018.formatted <- OH.2018.formatted[order(OH.2018.formatted$rs),]
#####################################################################################
# Section 3.0: Prepare to convert genotype files to numeric in GAPIT
# 3.1 Need to make column names first row in all three genotype files
Taxa <- colnames(C2.2018.formatted)
C2.2018.formatted <- rbind(Taxa,C2.2018.formatted)
Taxa <- colnames(OH.2018.formatted)
OH.2018.formatted <- rbind(Taxa, OH.2018.formatted)
Taxa <- colnames(SK.2016.formatted)
SK.2016.formatted <- rbind(Taxa, SK.2016.formatted)
# 3.2 Can remove first 12 columns because SNP names match SK dataframe!
C2.2018.formatted <-  dplyr::select(C2.2018.formatted, -(1:11))
# Need to add dplyr:: so R uses dplyr select rather than MASS
OH.2018.formatted <-  dplyr::select(OH.2018.formatted, -(1:11))
# 3.3 Merge all three genotype objects
Geno <- cbind(SK.2016.formatted,C2.2018.formatted,OH.2018.formatted)
####################################################################################
# Section 4.0: Convert to numeric format in GAPIT
# 4.1 Load GAPIT functions from source
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")  # 4.2 Create and move to directory $
dir.create("../GAPIT.Output/")
setwd("../GAPIT.Output/")
# 4.3 Run GAPIT numeric conversion
myGAPIT <- GAPIT(G=Geno, output.numerical=TRUE)
# writes numeric file to GAPIT.Output/ directory
#####################################################################################
# Section 5.0: Clean up R environment & Load Numeric Genotype
# 5.1 identify R objects no longer needed
Extras <- ls()
# 5.2 Remove those objects
rm(list=Extras)
rm(Extras)
# 5.3 read output files into R environment
Geno <- read_delim(file = "GAPIT.Genotype.Numerical.txt", col_names = TRUE, delim = "\t")
# 5.4 Move back to previous directory
setwd("../Rscripts")
#####################################################################################
# Section 6.0: Load phenotypic data to R environment
# 6.1 Read in phenotype data for all three experiments
C2.2018.Y <- read_delim(file="../Input.phenotype/C2.2018.pheno.csv",
                        col_names = T, delim = ",")
OH.2018.Y <- read_delim(file="../Input.phenotype/OH.2018.pheno.csv",
                        col_names = T, delim = ",")
SK.2016.Y <- read_delim(file="../Input.phenotype/SK.2016.pheno.csv",
                        col_names = T, delim = ",")
# 6.2 Combine all three dataframes
Pheno <- rbind(C2.2018.Y,OH.2018.Y,SK.2016.Y)
# 6.3 Match Phenotypic data to Genotypic data
Geno <- Geno[which(Geno$taxa %in% Pheno$taxa),]
Geno <- Geno[order(Geno$taxa),]
# 6.4 Match order of Phenotypic data to genotypic
Pheno <- Pheno[order(Pheno$taxa),]
# 6.4 Remove Extra R objects
rm(list=c("C2.2018.Y", "OH.2018.Y", "SK.2016.Y"))
#####################################################################################
# Section 7.0: Load Populations
# 7.1 Load population assignment file created during Rolling et al. 2020a
Exp.panel.IDs <- read_csv("../Input.data/Experiment_Panel.IDs.csv", col_names = TRUE)
# 7.2 Change "Taxa" to "taxa" so Exp.panel.IDs and Pheno can be merged
colnames(Exp.panel.IDs)[1] <- "taxa"
# 7.3 Load Population assignments into phenotype file via merge
Pheno_Popped <- merge(Pheno, Exp.panel.IDs, by="taxa")
# 7.4 Clean up R environment
rm(Exp.panel.IDs)
#####################################################################################
# 8.0 Define two options for Genomic Prediciton in BGLR
# 8.1 define bayesian iterations
nIter=60000
# 8.2 define burn in value
burnIn=20000
# 8.3 Make output matrix to store data
Output.data <- matrix(nrow=11, ncol=10)
# 8.4 Fill first column with Names of training and validation populations
Train.Valid.Names <- c("C2.US.GRIN", "PlusOne", "PlusTwo", "PlusThree",
                       "PlusFour", "PlusFive", "PlusSix","PlusSeven",
                       "US_like", "NotUSlike")
Output.data[1,] <- Train.Valid.Names
#####################################################################################
# Section 9.0 Confirm that C2.US predicts C2.GRIN very poorly!
# 9.1 Start a a loop of 10 iterations where C2.US predicts C2.GRIN
C2.US.Train.GRIN.valid <- lapply(1:10, function(x){
# 9.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 9.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 9.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
C2.GRIN.NA$ISW[1:153] <- NA
# 9.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 9.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 9.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 9.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 9.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 9.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
# 9.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[1:153])
return(Pearson.Cor$estimate) 
# 9.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 9.13 Close loop and function
  }
)
# 9.14 Store output in the Output.data object.
Output.data[2:11,1] <- unlist(C2.US.Train.GRIN.valid)
#####################################################################################
# Section 10.0 Add one GRIN sample to see if it improves things!
# 10.1 Start a loop to run 10 iterations
C2.US.TrainPlus1.GRIN.valid <- lapply(1:10, function(x){
# 10.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 10.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 10.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=1)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 10.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 10.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 10.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 10.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 10.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 10.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs]
# 10.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 10.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 10.13 Close loop and function
  }
)
# 10.14 Store output in the Output.data object.
Output.data[2:11,2] <- unlist(C2.US.TrainPlus1.GRIN.valid)
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
###############################################################
# Section 11.0. Add two GRIN sample to see if it improves things!
# 11.1 Start a loop to run 10 iterations
C2.US.TrainPlus2.GRIN.valid <- lapply(1:10, function(x){
# 11.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 11.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 11.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=2)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 11.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 11.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 11.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 11.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 11.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 11.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs]
# 11.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 11.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 11.13 Close loop and function
  }
)  
# 11.14 Store output in the Output.data object.
Output.data[2:11,3] <- unlist(C2.US.TrainPlus2.GRIN.valid)
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
##############################################################
# Section 12.0. Add three GRIN sample to see if it improves things!
# 12.1 Start a loop to run 10 iterations
C2.US.TrainPlus3.GRIN.valid <- lapply(1:10, function(x){
# 12.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 12.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 12.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=3)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 12.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 12.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 12.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 12.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 12.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 12.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs]
# 12.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 12.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 12.13 Close loop and function
  }
)  
# 12.14 Store output in the Output.data object.
Output.data[2:11,4] <- unlist(C2.US.TrainPlus3.GRIN.valid)
###############################################################
# Section 13.0. Add four GRIN sample to see if it improves things!
# 13.1 Start a loop to run 10 iterations
C2.US.TrainPlus4.GRIN.valid <- lapply(1:10, function(x){
# 13.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 13.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 13.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=4)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 13.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 13.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 13.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 13.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 13.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 13.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs]
# 13.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 13.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 13.13 Close loop and function
  }
)   
# 13.14 Store output in the Output.data object.
Output.data[2:11,5] <- unlist(C2.US.TrainPlus4.GRIN.valid)
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
###############################################################
# Section 14.0 Add five GRIN sample to see if it improves things!
# 14.1 Start a loop to run 10 iterations
C2.US.TrainPlus5.GRIN.valid <- lapply(1:10, function(x){
# 14.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 14.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 14.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=5)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 14.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 14.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 14.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 14.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 14.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 14.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs]
# 14.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 14.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 14.13 Close loop and function
  }
)  
# 14.14 Store output in the Output.data object.
Output.data[2:11,6] <- unlist(C2.US.TrainPlus5.GRIN.valid)
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
###############################################################
# Section 15.0. Add six GRIN sample to see if it improves things!
# 15.1 Start a loop to run 10 iterations
C2.US.TrainPlus6.GRIN.valid <- lapply(1:10, function(x){
# 15.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 15.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 15.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=6)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 15.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 15.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 15.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 15.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 15.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 15.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs] 
# 15.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 15.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 15.13 Close loop and function
  }
)  
# 15.14 Store output in the Output.data object.
Output.data[2:11,7] <- unlist(C2.US.TrainPlus6.GRIN.valid)
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
###############################################################
# Section 16.0 Add seven GRIN sample to see if it improves things! 
# 16.1 Start a loop to run 10 iterations
C2.US.TrainPlus7.GRIN.valid <- lapply(1:10, function(x){
# 16.2 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 16.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 16.4 Make object without C2.GRIN.Y ISW BLUEs
C2.GRIN.NA <- C2.GRIN.Y
random.sample <- sample(1:153, size=7)
frame <- 1:153
NAs <- setdiff(frame, random.sample)
C2.GRIN.NA$ISW[NAs] <- NA 
# 16.5 Merge dataframes
Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 16.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 16.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 16.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 16.9 Run BGLR Genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 16.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the $
dat <- dat[NAs]
# 16.11 Cor.test to compare predicted an observed values!
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW[NAs])
return(Pearson.Cor$estimate)
# 16.12 Clean up R environment
rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
# 16.13 Close loop and function
  }
)   
# 16.14 Store output in the Output.data object.
Output.data[2:11,8] <- unlist(C2.US.TrainPlus7.GRIN.valid)    
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
###############################################################
# Section 17.0 Add five GRIN samples with US identity... 
# 17.1 List Samples the select based on US population assignment
  # manually identified these PIs from Q-matrix
ListoUS <- c("PI592912A","PI548568","PI504486","PI592907A",
             "PI507686A")
# 17.2 Start a loop to run 10 iterations
C2.US.TrainPlusUS.GRIN.valid <- lapply(1:10, function(x){
# 17.3 Filter the Phenotype file to find PIs assigned to C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 17.4 Find those in the C2.GRIN panel
  C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# # 17.5 reorder dataframe so US-like PIs are at front  
  C2.GRIN_US.Y <- subset(C2.GRIN.Y, C2.GRIN.Y$taxa %in% ListoUS)
  C2.GRIN.2.Y <- subset(C2.GRIN.Y, !(C2.GRIN.Y$taxa %in% ListoUS))
  C2.GRIN2.Y <- rbind(C2.GRIN_US.Y, C2.GRIN.2.Y)
  C2.GRIN.NA <- C2.GRIN2.Y
# 17.6 Make object without C2.GRIN.Y ISW BLUEs
  C2.GRIN.NA$ISW[6:153] <- NA
# 17.7 Merge dataframes
  Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
# 17.8 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
# 17.9 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 17.10 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
# 17.11 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 17.12 Extract predicted phenotypes for C2.GRIN panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(C2.US.Y$taxa)+6):length(dat)] # subset to only consider those for the $
# 17.13 Cor.test to compare predicted an observed values!
  Pearson.Cor <- cor.test(dat, C2.GRIN2.Y$ISW[6:153])
  return(Pearson.Cor$estimate)
  # 17.14 Clean up R environment
  rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 17.15 Close loop and function
}
)   
# 17.16 Store output in the Output.data object.
Output.data[2:11,9] <- unlist(C2.US.TrainPlusUS.GRIN.valid)    

write.csv(x=Output.data, file="AddResUS2GRIN.csv")
###############################################################
# Section 18.0 Add five GRIN samples without US identity... 
# 18.1 List Samples the select based on US population assignment
# manually identified these PIs from Q-matrix
ListSelected <- c("PI189867","PI273483E","PI423823",
             "PI592912A","PI548568","PI504486")
# 18.2 Start a loop to run 10 iterations
C2.US.TrainPlusSelected.GRIN.valid <- lapply(1:10, function(x){
  # 18.3 Filter the Phenotype file to find PIs assigned to C2.US panel
  C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
  # 18.4 Find those in the C2.GRIN panel
  C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  # 18.5 reorder dataframe so US-like PIs are at front  
  C2.GRIN_GRIN.Y <- subset(C2.GRIN.Y, C2.GRIN.Y$taxa %in% ListSelected)
  C2.GRIN.2.Y <- subset(C2.GRIN.Y, !(C2.GRIN.Y$taxa %in% ListSelected))
  C2.GRIN.Y2 <- rbind(C2.GRIN_GRIN.Y, C2.GRIN.2.Y)
  C2.GRIN.NA <- C2.GRIN.Y2
  # 18.6 Make object without C2.GRIN.Y ISW BLUEs
  C2.GRIN.NA$ISW[7:153] <- NA
  # 18.7 Merge dataframes
  Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
  # 18.8 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 18.9 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 18.10 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
  # 18.11 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 18.12 Extract predicted phenotypes for C2.GRIN panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(C2.US.Y$taxa)+7):length(dat)] # subset to only consider those for the $
  # 18.13 Cor.test to compare predicted an observed values!
  Pearson.Cor <- cor.test(dat, C2.GRIN.Y2$ISW[7:153])
  return(Pearson.Cor$estimate)
  # 18.14 Clean up R environment
  rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 18.15 Close loop and function
}
)   
# 18.16 Store output in the Output.data object.
Output.data[2:11,10] <- unlist(C2.US.TrainPlusUS.GRIN.valid)    
write.csv(x=Output.data, file="AddResUS2GRIN.csv")

###############################################################
# Section 19.0 Add three US and three none US  
# 19.1 List Samples the select based on US population assignment
# manually identified these PIs from Q-matrix
ListGRIN <- c("PI189867","PI273483E","PI423823","PI458175A",
              "PI458244A")
# 19.2 Start a loop to run 10 iterations
C2.US.TrainPlusUS.GRIN.valid <- lapply(1:10, function(x){
  # 19.3 Filter the Phenotype file to find PIs assigned to C2.US panel
  C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
  # 19.4 Find those in the C2.GRIN panel
  C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  # 19.5 reorder dataframe so US-like PIs are at front  
  C2.GRIN_GRIN.Y <- subset(C2.GRIN.Y, C2.GRIN.Y$taxa %in% ListGRIN)
  C2.GRIN.2.Y <- subset(C2.GRIN.Y, !(C2.GRIN.Y$taxa %in% ListGRIN))
  C2.GRIN.Y2 <- rbind(C2.GRIN_GRIN.Y, C2.GRIN.2.Y)
  C2.GRIN.NA <- C2.GRIN.Y2
  # 19.6 Make object without C2.GRIN.Y ISW BLUEs
  C2.GRIN.NA$ISW[6:153] <- NA
  # 19.7 Merge dataframes
  Exp.Y <- rbind(C2.US.Y, C2.GRIN.NA)
  # 19.8 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 19.9 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 19.10 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]
  # 19.11 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 19.12 Extract predicted phenotypes for C2.GRIN panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(C2.US.Y$taxa)+6):length(dat)] # subset to only consider those for the $
  # 19.13 Cor.test to compare predicted an observed values!
  Pearson.Cor <- cor.test(dat, C2.GRIN.Y2$ISW[6:153])
  return(Pearson.Cor$estimate)
  # 19.14 Clean up R environment
  rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 19.15 Close loop and function
}
)   
# 19.16 Store output in the Output.data object.
Output.data[2:11,10] <- unlist(C2.US.TrainPlusUS.GRIN.valid) 
write.csv(x=Output.data, file="AddResUS2GRIN.csv")
# The End :) 
quit()
