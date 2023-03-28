# Results:Genomic prediction using SK PI panels to generate or validate models
# Author: William Rolling
# Corresponding Author: Leah McHale (mchale.21@osu.edu)
#####################################################################################
# Section 1.0: Prepare R environment
# 1.1 Clear R environment
rm(list=ls()) 
# 1.2 List packages needed to complete analysis
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
# example for someone who downloads data to folder in desktop <---- Set pathway to folder here. 
# setwd("~/../Desktop/Supplemental Data 1/")
#####################################################################################
# Section 2.0: Load genotypic data into R environment    
# 2.1 Load genotype files for three experiments
SK.2016 <- read_csv("Input.data/SK.2016.hmp.csv", col_names=TRUE)
C2.2018 <- read_csv("Input.data/C2.2018.hmp.csv", col_names=TRUE)
OH.2018 <- read_csv("Input.data/OH.2018.hmp.csv", col_names=TRUE)
# 2.2 rename first column to match SK.2016 column name
colnames(C2.2018)[1] <- "rs"
colnames(OH.2018)[1] <- "rs"
colnames(SK.2016)[1] <- "rs"

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
#####################################################################################
# Section 4.0: Convert to numeric format in GAPIT     
# 4.1 Load GAPIT functions from source  
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")  # 4.2 Create and move to directory to store the following output
dir.create("GAPIT.Output/")
setwd("GAPIT.Output/")
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
  setwd("../")
#####################################################################################
# Section 6.0: Load phenotypic data to R environment 
# 6.1 Read in phenotype data for all three experiments  
  C2.2018.Y <- read_delim(file="Input.phenotype/C2.2018.pheno.csv", col_names = T, delim = ",") 
  OH.2018.Y <- read_delim(file="Input.phenotype/OH.2018.pheno.csv", col_names = T, delim = ",")   
  SK.2016.Y <- read_delim(file="Input.phenotype/SK.2016.pheno.csv", col_names = T, delim = ",")  
  # 6.2 Combine all three dataframes
  Pheno <- rbind(C2.2018.Y,OH.2018.Y,SK.2016.Y)
# 6.3 Match Genotypic data to phenotypic data
Geno <- Geno[which(Geno$taxa %in% Pheno$taxa),]
Geno <- Geno[order(Geno$taxa),]
# 6.4 Match Phenotypic data to genotypic
Pheno <- Pheno[order(Pheno$taxa),]
# 6.5 Remove Extra R objects
rm(list=c("C2.2018.Y", "OH.2018.Y", "SK.2016.Y"))
#####################################################################################
# Section 7.0: Load Populations
# 7.1 Load population assignment file created during Rolling et al. 2020a
Exp.panel.IDs <- read_csv("Input.data/Experiment_Panel.IDs.csv", col_names = TRUE)
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
Output.data <- matrix(nrow=8, ncol=4)
# 8.4 Fill first column with Names of training and validation populations  
Train.Valid.Names <- c("2016-SK Train and OH.SK Valid",
                       "2016-SK Train and C2.SK Valid",
                       "2016-SK Train and OH.US Valid",
                       "2016-SK Train and C2.US Valid",
                       "2016-SK Train and OH.GRIN Valid",
                       "2016-SK Train and C2.GRIN Valid",
                       "2016-SK Train and 2018-OH Valid",
                       "2016-SK Train and 2018-C2 Valid")
Output.data[1:8,1] <- Train.Valid.Names 
#####################################################################################      
# Section 9.0: R.SK as training population for OH.SK panel
# 9.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
 R.SKTrain.OH.SK.Valid <- lapply(1:10, function(x){
  # 9.2 Filter the Phenotype file to find PIs assigned to OH.SK panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
  # 9.3 Find those in the OH.SK panel
  OH.SK.Y <- filter(Pheno_Popped, Population == "OH.SK")
  # 9.4 Make object without OH.SK ISW BLUEs
  OH.SK.Y.NA <- OH.SK.Y
  OH.SK.Y.NA$ISW <- NA
  # 9.5 Merge dataframes  
  Exp.Y <- rbind(R.SK.Y, OH.SK.Y.NA)
  # 9.6 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 9.7 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 9.8 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 9.9 Run BGLR genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 9.10 Extract predicted phenotypes for OH.SK panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel
  # 9.11 Cor.test to compare predicted an observed values! 
  Pearson.Cor <- cor.test(dat, OH.SK.Y$ISW)
  return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
  # 9.12 Clean up R environment
  rm(list= c("OH.SK.Y", "OH.SK.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 9.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.OH.SK.Valid <- as.numeric(unlist(R.SKTrain.OH.SK.Valid))
# 9.14  Add average to output data frame and the standard deviation
Output.data[1,2] <- mean(R.SKTrain.OH.SK.Valid)
Output.data[1,3] <- sd(R.SKTrain.OH.SK.Valid)

write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 10.0: R.SK as training population for C2.SK panel
# 10.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.C2.SK.Valid <- lapply(1:10, function(x){
  # 10.2 Filter the Phenotype file to find PIs assigned to C2.SK panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 10.3 Find those in the C2.SK panel
C2.SK.Y <- filter(Pheno_Popped, Population == "C2.SK")
# 10.4 Make object without C2.SK ISW BLUEs
C2.SK.Y.NA <- C2.SK.Y
C2.SK.Y.NA$ISW <- NA
# 10.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, C2.SK.Y.NA)
# 10.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 10.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 10.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 10.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 10.10 Extract predicted phenotypes for C2.SK panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.SK panel
# 10.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, C2.SK.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 10.12 Clean up R environment
rm(list= c("C2.SK.Y", "C2.SK.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 10.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.C2.SK.Valid <- as.numeric(unlist(R.SKTrain.C2.SK.Valid))
# 10.14  Add average to output data frame and the standard deviation
Output.data[2,2] <- mean(R.SKTrain.C2.SK.Valid)
Output.data[2,3] <- sd(R.SKTrain.C2.SK.Valid)
write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 11.0: R.SK as training population for OH.US panel
# 11.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.OH.US.Valid <- lapply(1:10, function(x){
  # 11.2 Filter the Phenotype file to find PIs assigned to OH.US panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 11.3 Find those in the OH.US pane
OH.US.Y <- filter(Pheno_Popped, Population == "OH.US")
# 11.4 Make object without OH.US ISW BLUEs
OH.US.Y.NA <- OH.US.Y
OH.US.Y.NA$ISW <- NA
# 11.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, OH.US.Y.NA)
# 11.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 11.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 11.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 11.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 11.10 Extract predicted phenotypes for OH.US panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.US panel
# 11.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, OH.US.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 11.12 Clean up R environment
rm(list= c("OH.US.Y", "OH.US.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 11.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.OH.US.Valid <- as.numeric(unlist(R.SKTrain.OH.US.Valid))
# 11.14  Add average to output data frame and the standard deviation
Output.data[3,2] <- mean(R.SKTrain.OH.US.Valid)
Output.data[3,3] <- sd(R.SKTrain.OH.US.Valid)
write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 12.0: R.SK as training population for C2.US panel
# 12.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.C2.US.Valid <- lapply(1:10, function(x){
  # 12.2 Filter the Phenotype file to find PIs assigned to C2.US panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 12.3 Find those in the C2.US panel
C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
# 12.4 Make object without C2.US ISW BLUEs
C2.US.Y.NA <- C2.US.Y
C2.US.Y.NA$ISW <- NA
# 12.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, C2.US.Y.NA)
# 12.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 12.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 12.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 12.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 12.10 Extract predicted phenotypes for C2.US panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.US panel
# 12.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, C2.US.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 12.12 Clean up R environment
rm(list= c("C2.US.Y", "C2.US.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 12.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.C2.US.Valid <- as.numeric(unlist(R.SKTrain.C2.US.Valid))
# 12.14  Add average to output data frame and the standard deviation
Output.data[4,2] <- mean(R.SKTrain.C2.US.Valid)
Output.data[4,3] <- sd(R.SKTrain.C2.US.Valid)
write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 13.0: R.SK as training population for OH.GRIN panel
# 13.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.OH.GRIN.Valid <- lapply(1:10, function(x){
  # 13.2 Filter the Phenotype file to find PIs assigned to OH.GRIN panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 13.3 Find those in the OH.GRIN panel
OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")
# 13.4 Make object without OH.GRIN ISW BLUEs
OH.GRIN.Y.NA <- OH.GRIN.Y
OH.GRIN.Y.NA$ISW <- NA
# 13.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, OH.GRIN.Y.NA)
# 13.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 13.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 13.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 13.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 13.10 Extract predicted phenotypes for OH.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.GRIN panel
# 13.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, OH.GRIN.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 13.12 Clean up R environment
rm(list= c("OH.GRIN.Y", "OH.GRIN.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 13.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.OH.GRIN.Valid <- as.numeric(unlist(R.SKTrain.OH.GRIN.Valid))
# 13.14  Add average to output data frame and the standard deviation
Output.data[5,2] <- mean(R.SKTrain.OH.GRIN.Valid)
Output.data[5,3] <- sd(R.SKTrain.OH.GRIN.Valid)
write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 14.0: R.SK as training population for C2.GRIN panel
# 14.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.C2.GRIN.Valid <- lapply(1:10, function(x){
  # 14.2 Filter the Phenotype file to find PIs assigned to C2.GRIN panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 14.3 Find those in the C2.GRIN panel
C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
# 14.4 Make object without C2.GRIN ISW BLUEs
C2.GRIN.Y.NA <- C2.GRIN.Y
C2.GRIN.Y.NA$ISW <- NA
# 14.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, C2.GRIN.Y.NA)
# 14.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 14.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 14.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 14.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 14.10 Extract predicted phenotypes for C2.GRIN panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.GRIN panel
# 14.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 14.12 Clean up R environment
rm(list= c("C2.GRIN.Y", "C2.GRIN.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 14.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.C2.GRIN.Valid <- as.numeric(unlist(R.SKTrain.C2.GRIN.Valid))
# 14.14  Add average to output data frame and the standard deviation
Output.data[6,2] <- mean(R.SKTrain.C2.GRIN.Valid)
Output.data[6,3] <- sd(R.SKTrain.C2.GRIN.Valid)
write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 15.0: R.SK as training population for OH2018 panel
# 15.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.OH2018.Valid <- lapply(1:10, function(x){
  # 15.2 Filter the Phenotype file to find PIs assigned to OH2018 panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 15.3 Find those in the OH2018 panel
OH2018.Y <- filter(Pheno_Popped, Population %in% c("OH.SK", "OH.US", "OH.GRIN"))
# 15.4 Make object without OH2018 ISW BLUEs
OH2018.Y.NA <- OH2018.Y
OH2018.Y.NA$ISW <- NA
# 15.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, OH2018.Y.NA)
# 15.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 15.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 15.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 15.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 15.10 Extract predicted phenotypes for OH2018 panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the OH2018 panel
# 15.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, OH2018.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 15.12 Clean up R environment
rm(list= c("OH2018.Y", "OH2018.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 15.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.OH2018.Valid <- as.numeric(unlist(R.SKTrain.OH2018.Valid))
# 15.14  Add average to output data frame and the standard deviation
Output.data[7,2] <- mean(R.SKTrain.OH2018.Valid)
Output.data[7,3] <- sd(R.SKTrain.OH2018.Valid)
write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)
#####################################################################################      
# Section 16.0: R.SK as training population for C2018 panel
# 16.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average 
R.SKTrain.C2018.Valid <- lapply(1:10, function(x){
  # 16.2 Filter the Phenotype file to find PIs assigned to C2018 panel
  R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
# 16.3 Find those in the C2018 panel
C2018.Y <- filter(Pheno_Popped, Population %in% c("C2.SK", "C2.US", "C2.GRIN"))
# 16.4 Make object without C2018 ISW BLUEs
C2018.Y.NA <- C2018.Y
C2018.Y.NA$ISW <- NA
# 16.5 Merge dataframes  
Exp.Y <- rbind(R.SK.Y, C2018.Y.NA)
# 16.6 Remove Population Assignment
Exp.Y <- Exp.Y[,-3]
# 16.7 Make genotype file
Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
# 16.8 Match order between dataframes!
Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
# 16.9 Run BGLR genomic Prediction
ETA<-list(list(X=Exp.G[,-1], model='BL'))
fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
# 16.10 Extract predicted phenotypes for C2018 panel
dat <- fm$yHat # extracts predicted phenotypes from BGLR output
dat <- dat[(length(R.SK.Y$taxa)+1):length(dat)] # subset to only consider those for the C2018 panel
# 16.11 Cor.test to compare predicted an observed values! 
Pearson.Cor <- cor.test(dat, C2018.Y$ISW)
return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
# 16.12 Clean up R environment
rm(list= c("C2018.Y", "C2018.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
}
)

# 16.13 unlist and convert to numeric so I can calculate average  
R.SKTrain.C2018.Valid <- as.numeric(unlist(R.SKTrain.C2018.Valid))
# 16.14  Add average to output data frame and the standard deviation
Output.data[8,2] <- mean(R.SKTrain.C2018.Valid)
Output.data[8,3] <- sd(R.SKTrain.C2018.Valid)

write.csv(x=Output.data, file="SK2016Trainer.csv", row.names=FALSE)