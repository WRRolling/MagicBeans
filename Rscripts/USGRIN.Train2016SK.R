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
# Create dataframe for US & GRIN panels to train 2016 SK panels
# 8.3 Make output matrix to store data 
Output.data <- matrix(nrow=5, ncol=3)
# 8.4 Fill first column with Names of training and validation populations  
Train.Valid.Names <- c("C2.GRIN Train and 2016-SK Valid",
                       "C2.US Train and 2016-SK Valid",
                       "OH.GRIN and 2016-SK Valid",
                       "OH.US Train and 2016-SK Valid") 
# 8.5 Add names to output data
Output.data[1:4,1] <- Train.Valid.Names
#####################################################################################          
# Section 9.0: C2.GRIN as training population for 2016.SK panel  
# 9.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average   
C2.GRINTrain.2016.SK.Valid <- lapply(1:10, function(x){
  # 9.2 Find those in the OH.SK panel
  C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  C2.2016.SK.Y <- filter(Pheno_Popped, Population == "R.SK") 
  # 9.3 Make object without ISW BLUEs
  C2.2016.SK.Y.NA <- C2.2016.SK.Y
  C2.2016.SK.Y.NA$ISW <- NA
  # 9.4 Merge dataframes  
  Exp.Y <- rbind(C2.GRIN.Y, C2.2016.SK.Y.NA)
  # 9.5 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 9.6 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 9.7 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 9.8 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 9.9 Extract predicted phenotypes for OH.SK panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(C2.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel  
  # 9.10 Cor.test to compare predicted an observed values! 
  Pearson.Cor <- cor.test(dat, C2.2016.SK.Y$ISW)
  return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 9.11 Clean up R environment
  rm(list= c("C2.2016.SK.Y", "C2.GRIN.Y", "C2.2016.SK.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))    
}
)
# 9.12 unlist and convert to numeric so I can calculate average  
C2.GRINTrain.2016.SK.Valid  <- as.numeric(unlist(C2.GRINTrain.2016.SK.Valid ))
# 9.13 Add average to output dataframe and the standard deviation
Output.data[1,2] <- mean(C2.GRINTrain.2016.SK.Valid )
Output.data[1,3] <- sd(C2.GRINTrain.2016.SK.Valid )
#####################################################################################          
# Section 10.0: C2.US as training population for 2016.SK panel  
# 10.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average       
C2.USTrain.2016.SK.Valid <- lapply(1:10, function(x){
  # 10.2 Find those in the OH.SK panel
  C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
  C2.2016.SK.Y <- filter(Pheno_Popped, Population == "R.SK") 
  # 10.3 Make object without ISW BLUEs
  C2.2016.SK.Y.NA <- C2.2016.SK.Y
  C2.2016.SK.Y.NA$ISW <- NA
  # 10.4 Merge dataframes  
  Exp.Y <- rbind(C2.US.Y, C2.2016.SK.Y.NA)
  # 10.5 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 10.6 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 10.7 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 10.8 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 10.9 Extract predicted phenotypes for OH.SK panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(C2.US.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel  
  # 10.10 Cor.test to compare predicted an observed values! 
  Pearson.Cor <- cor.test(dat, C2.2016.SK.Y$ISW)
  return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 10.11 Clean up R environment
  rm(list= c("C2.2016.SK.Y", "C2.US.Y", "C2.2016.SK.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))    
}
)
# 10.12 unlist and convert to numeric so I can calculate average  
C2.USTrain.2016.SK.Valid  <- as.numeric(unlist(C2.USTrain.2016.SK.Valid ))
# 10.13 Add average to output dataframe and the standard deviation
Output.data[2,2] <- mean(C2.USTrain.2016.SK.Valid )
Output.data[2,3] <- sd(C2.USTrain.2016.SK.Valid )    
####################################################################################          
# Section 11.0: OH.GRIN as training population for 2016.SK panel  
# 11.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average       
OH.GRINTrain.2016.SK.Valid <- lapply(1:10, function(x){
  # 11.2 Find those in the OH.SK panel
  OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")
  SK.2016.Y <- filter(Pheno_Popped, Population == "R.SK") 
  # 11.3 Make object without ISW BLUEs
  SK.2016.Y.NA <- SK.2016.Y
  SK.2016.Y.NA$ISW <- NA
  # 11.4 Merge dataframes  
  Exp.Y <- rbind(OH.GRIN.Y, SK.2016.Y.NA)
  # 11.5 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 11.6 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 11.7 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 11.8 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 11.9 Extract predicted phenotypes for OH.SK panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(OH.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel  
  # 11.10 Cor.test to compare predicted an observed values! 
  Pearson.Cor <- cor.test(dat, SK.2016.Y$ISW)
  return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 11.11 Clean up R environment
  rm(list= c("OH.GRIN.Y", "SK.2016.Y", "SK.2016.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))    
}
)
# 11.12 unlist and convert to numeric so I can calculate average  
OH.GRINTrain.2016.SK.Valid  <- as.numeric(unlist(OH.GRINTrain.2016.SK.Valid ))
# 11.13 Add average to output dataframe and the standard deviation
Output.data[3,2] <- mean(OH.GRINTrain.2016.SK.Valid)
Output.data[3,3] <- sd(OH.GRINTrain.2016.SK.Valid)            

####################################################################################          
# Section 12.0: OH.GRIN as training population for 2016.SK panel  
# 12.1 Start a loop. 
# The Gibbs sampler will result in slightly different estimations.
# As such I will run 10 iterations and get an average       
OH.USTrain.2016.SK.Valid <- lapply(1:10, function(x){
  # 12.2 Find those in the OH.SK panel
  OH.US.Y <- filter(Pheno_Popped, Population == "OH.US")
  SK.2016.Y <- filter(Pheno_Popped, Population == "R.SK") 
  # 12.3 Make object without ISW BLUEs
  SK.2016.Y.NA <- SK.2016.Y
  SK.2016.Y.NA$ISW <- NA
  # 12.4 Merge dataframes  
  Exp.Y <- rbind(OH.US.Y, SK.2016.Y.NA)
  # 12.5 Remove Population Assignment
  Exp.Y <- Exp.Y[,-3]
  # 12.6 Make genotype file
  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 12.7 Match order between dataframes!
  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 12.8 Run BGLR Genomic Prediction
  ETA<-list(list(X=Exp.G[,-1], model='BL'))
  fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 12.9 Extract predicted phenotypes for OH.SK panel
  dat <- fm$yHat # extracts predicted phenotypes from BGLR output
  dat <- dat[(length(OH.US.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel  
  # 12.10 Cor.test to compare predicted an observed values! 
  Pearson.Cor <- cor.test(dat, SK.2016.Y$ISW)
  return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 12.11 Clean up R environment
  rm(list= c("OH.US.Y", "SK.2016.Y", "SK.2016.Y.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))    
}
)
# 12.12 unlist and convert to numeric so I can calculate average  
OH.USTrain.2016.SK.Valid  <- as.numeric(unlist(OH.USTrain.2016.SK.Valid ))
# 12.13 Add average to output dataframe and the standard deviation
Output.data[4,2] <- mean(OH.USTrain.2016.SK.Valid)
Output.data[4,3] <- sd(OH.USTrain.2016.SK.Valid)            

#####################################################################################      
# Section Thirteen: Write Double Training Panel Output to file 
write.csv(x=Output.data, file="US.GRIN.Train2016SK.csv", row.names = F)

#####################################################################################  
