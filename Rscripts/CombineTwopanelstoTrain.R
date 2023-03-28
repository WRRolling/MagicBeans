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
    Output.data <- matrix(nrow=3, ncol=3)
  # 8.4 Fill first column with Names of training and validation populations  
    Train.Valid.Names <- c("C2.SK.US Train and C2.GRIN Valid",
                           "C2.SK.GRIN Train and C2.US Valid",
                           "C2.US.GRIN Train and C2.SK Valid")
    Output.data[1:3,1] <- Train.Valid.Names
    
#####################################################################################      
# Section 9.0: C2.SK & US panels as training population for C2.GRIN panel
  # 9.1 Start a loop. 
    # The Gibbs sampler will result in slightly different estimations.
    # As such I will run 10 iterations and get an average 
    C2.SKUSTrain.GRIN.Valid <- lapply(1:10, function(x){
  # 9.2 Filter the Phenotype file to find PIs assigned to C2.SK panel
      C2.SK.US.Y <- filter(Pheno_Popped, Population %in% c("C2.SK","C2.US"))
  # 9.3 Find those in the C2.US panel
        C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  # 9.4 Make object without C2.US.Y ISW BLUEs
          C2.GRIN.NA <- C2.GRIN.Y
            C2.GRIN.NA$ISW <- NA
  # 9.5 Merge dataframes  
              Exp.Y <- rbind(C2.SK.US.Y, C2.GRIN.NA)
  # 9.6 Remove Population Assignment
                Exp.Y <- Exp.Y[,-3]
  # 9.7 Make genotype file
                  Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 9.8 Match order between dataframes!
                    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 9.9 Run BGLR genomic Prediction
                  ETA<-list(list(X=Exp.G[,-1], model='BL'))
                fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 9.10 Extract predicted phenotypes for C2.US panel
              dat <- fm$yHat # extracts predicted phenotypes from BGLR output
            dat <- dat[(length(C2.SK.US.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.US panel
  # 9.11 Cor.test to compare predicted an observed values! 
          Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW)
          return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
  # 9.12 Clean up R environment
      rm(list= c("C2.SK.US.Y", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
    }
    )
    
  # 9.13 unlist and convert to numeric so I can calculate average  
    C2.SKUSTrain.GRIN.Valid <- as.numeric(unlist(C2.SKUSTrain.GRIN.Valid))
  # 9.14  Add average to output data frame and the standard deviation
    Output.data[1,2] <- mean(C2.SKUSTrain.GRIN.Valid)
    Output.data[1,3] <- sd(C2.SKUSTrain.GRIN.Valid)

    
    
#####################################################################################      
# Section 10.0: C2.SK and GRIN panel as training population for C2.US panel  
    # 10.1 Start a loop. 
    # The Gibbs sampler will result in slightly different estimations.
    # As such I will run 10 iterations and get an average   
    C2.SK.GRIN.Train.US.Valid <- lapply(1:10, function(x){
  # 10.2 Find those in the C2.GRIN panel
    C2.SK.GRIN.Y <- filter(Pheno_Popped, Population %in% c("C2.SK", "C2.GRIN"))
    C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
  # 10.3 Make object without C2.GRIN.Y ISW BLUEs
    C2.US.NA <- C2.US.Y
    C2.US.NA$ISW <- NA
  # 10.4 Merge dataframes  
    Exp.Y <- rbind(C2.SK.GRIN.Y, C2.US.NA)
  # 10.5 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 10.6 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 10.7 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 10.8 Run BGLR Genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 10.9 Extract predicted phenotypes for C2.GRIN panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(C2.SK.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.GRIN panel
  # 10.10 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, C2.US.Y$ISW)
    return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 10.11 Clean up R environment
    rm(list= c("C2.SK.GRIN.Y", "C2.US.Y", "C2.US.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))   
      }
    )
  # 10.12 unlist and convert to numeric so I can calculate average  
    C2.SK.GRIN.Train.US.Valid <- as.numeric(unlist(C2.SK.GRIN.Train.US.Valid))
  # 10.13 Add average to output dataframe and the standard deviation
    Output.data[2,2] <- mean(C2.SK.GRIN.Train.US.Valid)
    Output.data[2,3] <- sd(C2.SK.GRIN.Train.US.Valid)

    
#####################################################################################      
# Section 11.0: C2.GRIN & US panels as training population for C2.SK panel  
  # 10.1 Start a loop. 
    # The Gibbs sampler will result in slightly different estimations.
    # As such I will run 10 iterations and get an average   
    C2.US.GRINTrain.SK.Valid <- lapply(1:10, function(x){
  # 11.2 Find those in the OH.SK panel
    C2.US.GRIN.Y <- filter(Pheno_Popped, Population %in% c("C2.US", "C2.GRIN"))
    C2.SK.Y <- filter(Pheno_Popped, Population == "C2.SK")
  # 11.3 Make object without OH.SK.Y ISW BLUEs
    C2.SK.NA <- C2.SK.Y
    C2.SK.NA$ISW <- NA
  # 11.4 Merge dataframes  
    Exp.Y <- rbind(C2.US.GRIN.Y, C2.SK.NA)
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
    dat <- dat[(length(C2.US.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel
  # 11.10 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, C2.SK.Y$ISW)
    return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 11.11 Clean up R environment
    rm(list= c("OH.SK.NA", "OH.SK.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))    
    }
    )
    # 11.12 unlist and convert to numeric so I can calculate average  
    C2.US.GRINTrain.SK.Valid  <- as.numeric(unlist(C2.US.GRINTrain.SK.Valid ))
    # 11.13 Add average to output data frame and the standard deviation
    Output.data[3,2] <- mean(C2.US.GRINTrain.SK.Valid )
    Output.data[3,3] <- sd(C2.US.GRINTrain.SK.Valid )

#####################################################################################      
# Section Twelve: Write Double Training Panel Output to file 
    write.csv(x=Output.data, file="C2Double.Train.csv", row.names = F)
    
#####################################################################################  
# 8.0 Define two options for Genomic Prediction in BGLR 
# 8.1 Make New output matrix to store data 
    Output.data <- matrix(nrow=3, ncol=4)
    # 8.4 Fill first column with Names of training and validation populations  
    Train.Valid.Names <- c("OH.SK.US Train and OH.GRIN Valid",
                           "OH.SK.GRIN Train and OH.US Valid",
                           "OH.US.GRIN Train and OH.SK Valid")
    Output.data[1:3,1] <- Train.Valid.Names
    
    #####################################################################################      
    # Section 9.0: OH.SK & US panels as training population for OH.GRIN panel
    # 9.1 Start a loop. 
    # The Gibbs sampler will result in slightly different estimations.
    # As such I will run 10 iterations and get an average 
    OH.SKUSTrain.GRIN.Valid <- lapply(1:10, function(x){
      # 9.2 Filter the Phenotype file to find PIs assigned to OH.SK and US panel
      OH.SK.US.Y <- filter(Pheno_Popped, Population %in% c("OH.SK","OH.US"))
      # 9.3 Find those in the OH.GRIN panels
      OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")
      # 9.4 Make object without OH.GRIN.Y ISW BLUEs
      OH.GRIN.NA <- OH.GRIN.Y
      OH.GRIN.NA$ISW <- NA
      # 9.5 Merge dataframes  
      Exp.Y <- rbind(OH.SK.US.Y, OH.GRIN.NA)
      # 9.6 Remove Population Assignment
      Exp.Y <- Exp.Y[,-3]
      # 9.7 Make genotype file
      Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
      # 9.8 Match order between dataframes!
      Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
      # 9.9 Run BGLR genomic Prediction
      ETA<-list(list(X=Exp.G[,-1], model='BL'))
      fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
      # 9.10 Extract predicted phenotypes for OH.GRI panel
      dat <- fm$yHat # extracts predicted phenotypes from BGLR output
      dat <- dat[(length(OH.SK.US.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.GRIN panel
      # 9.11 Cor.test to compare predicted an observed values! 
      Pearson.Cor <- cor.test(dat, OH.GRIN.Y$ISW)
      return(Pearson.Cor$estimate) # 0.70 average heritability of ISW
      # 9.12 Clean up R environment
      rm(list= c("OH.SK.US.Y", "OH.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
    }
    )
    
    # 9.13 unlist and convert to numeric so I can calculate average  
    OH.SKUSTrain.GRIN.Valid <- as.numeric(unlist(OH.SKUSTrain.GRIN.Valid))
    # 9.14  Add average to output data frame and the standard deviation
    Output.data[1,2] <- mean(OH.SKUSTrain.GRIN.Valid)
    Output.data[1,3] <- sd(OH.SKUSTrain.GRIN.Valid)
    
    
    
#####################################################################################      
    # Section 10.0: OH.SK and OH.GRIN panel as training population for OH.US panel  
    # 10.1 Start a loop. 
    # The Gibbs sampler will result in slightly different estimations.
    # As such I will run 10 iterations and get an average   
    OH.SK.GRIN.Train.US.Valid <- lapply(1:10, function(x){
      # 10.2 Find those in the SK and GRIN panel
      OH.SK.GRIN.Y <- filter(Pheno_Popped, Population %in% c("OH.SK", "OH.GRIN"))
      OH.US.Y <- filter(Pheno_Popped, Population == "OH.US")
      # 10.3 Make object without OH.US ISW BLUEs
      OH.US.NA <- OH.US.Y
      OH.US.NA$ISW <- NA
      # 10.4 Merge dataframes  
      Exp.Y <- rbind(OH.SK.GRIN.Y, OH.US.NA)
      # 10.5 Remove Population Assignment
      Exp.Y <- Exp.Y[,-3]
      # 10.6 Make genotype file
      Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
      # 10.7 Match order between dataframes!
      Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
      # 10.8 Run BGLR Genomic Prediction
      ETA<-list(list(X=Exp.G[,-1], model='BL'))
      fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
      # 10.9 Extract predicted phenotypes for OH.US panel
      dat <- fm$yHat # extracts predicted phenotypes from BGLR output
      dat <- dat[(length(OH.SK.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.US panel
      # 10.10 Cor.test to compare predicted an observed values! 
      Pearson.Cor <- cor.test(dat, OH.US.Y$ISW)
      return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
      # 10.11 Clean up R environment
      rm(list= c("OH.SK.GRIN.Y", "OH.US.Y", "OH.US.NA", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))   
    }
    )
    # 10.12 unlist and convert to numeric so I can calculate average  
    OH.SK.GRIN.Train.US.Valid <- as.numeric(unlist(OH.SK.GRIN.Train.US.Valid))
    # 10.13 Add average to output dataframe and the standard deviation
    Output.data[2,2] <- mean(OH.SK.GRIN.Train.US.Valid)
    Output.data[2,3] <- sd(OH.SK.GRIN.Train.US.Valid)
    
    
    #####################################################################################      
    # Section 11.0: OH.GRIN & US panels as training population for OH.SK panel  
    # 10.1 Start a loop. 
    # The Gibbs sampler will result in slightly different estimations.
    # As such I will run 10 iterations and get an average   
    OH.US.GRINTrain.SK.Valid <- lapply(1:10, function(x){
      # 11.2 Find those in the US and GRIN panels
      OH.US.GRIN.Y <- filter(Pheno_Popped, Population %in% c("OH.US", "OH.GRIN"))
      OH.SK.Y <- filter(Pheno_Popped, Population == "OH.SK")
      # 11.3 Make object without OH.SK.Y ISW BLUEs
      OH.SK.NA <- OH.SK.Y
      OH.SK.NA$ISW <- NA
      # 11.4 Merge dataframes  
      Exp.Y <- rbind(OH.US.GRIN.Y, OH.SK.NA)
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
      dat <- dat[(length(OH.US.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel
      # 11.10 Cor.test to compare predicted an observed values! 
      Pearson.Cor <- cor.test(dat, OH.SK.Y$ISW)
      return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
      # 11.11 Clean up R environment
      rm(list= c("OH.SK.NA", "OH.SK.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))    
    }
    )
    # 11.12 unlist and convert to numeric so I can calculate average  
    OH.US.GRINTrain.SK.Valid  <- as.numeric(unlist(C2.US.GRINTrain.SK.Valid ))
    # 11.13 Add average to output dataframe and the standard deviation
    Output.data[3,2] <- mean(OH.US.GRINTrain.SK.Valid )
    Output.data[3,3] <- sd(OH.US.GRINTrain.SK.Valid )
    
#####################################################################################      
    # Section Twelve: Write Double Training Panel Output to file 
    write.csv(x=Output.data, file="OHDouble.Train.csv", row.names = F)
#####################################################################################  
    
    
    
        