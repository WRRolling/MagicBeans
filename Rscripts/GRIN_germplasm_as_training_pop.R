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
    # MASS and dplyr 'select' do not play well together
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
    setwd("~/../Desktop/Supplemental Data 1/")
#####################################################################################
# Section 2.0: Load genotypic data into R environment    
  # 2.1 Load genotype files for three experiments
    SK.2016 <- read_csv("Input.data/SK.2016.hmp.csv", col_names=TRUE)
    C2.2018 <- read_csv("Input.data/C2.2018.hmp.csv", col_names=TRUE)
    OH.2018 <- read_csv("Input.data/OH.2018.hmp.csv", col_names=TRUE)
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
    Pheno <-  Pheno[-which(duplicated(Pheno$taxa)=="TRUE"),]
  # 6.3 Match Phenotypic data to Genotypic data
    Geno <- Geno[which(Geno$taxa %in% Pheno$taxa),]
    Geno <- Geno[order(Geno$taxa),]
  # 6.4 Match order of Phenotypic data to genotypic
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
    Output.data <- matrix(nrow=4, ncol=3)
  # 8.4 Fill first column with Names of training and validation populations  
    Train.Valid.Names <- c("C2.GRIN Train and C2.US Valid",
                           "C2.GRIN Train and C2.SK Valid",
                           "C2.GRIN Train and OH.GRIN Valid",
                           "C2.GRIN LOOCV")
    Output.data[1:4,1] <- Train.Valid.Names  
#####################################################################################      
# Section 9.0: C2.GRIN panel as training population for C2.US panel
  # 9.1 Start a loop so ten iterations can be completed
    C2.GRIN.Train.US.valid <- lapply(1:10, function(x){
  # 9.2 Filter the Phenotype file to find PIs assigned to C2.GRIN panel
    C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  # 9.3 Filter the Phenotype file to find PIs assigned to C2.US panel
      C2.US.Y <- filter(Pheno_Popped, Population == "C2.US")
  # 9.4 Make object without C2.US.Y ISW BLUEs
        C2.US.NA <- C2.US.Y
          C2.US.NA$ISW <- NA
  # 9.5 Merge dataframes  
            Exp.Y <- rbind(C2.GRIN.Y, C2.US.NA)
  # 9.6 Remove Population Assignment
              Exp.Y <- Exp.Y[,-3]
  # 9.7 Prepare genotype file
                Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 9.8 Match order between dataframes!
                  Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 9.9 Run BGLR genomic Prediction
                  ETA<-list(list(X=Exp.G[,-1], model='BL'))
                fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 9.10 Extract predicted phenotypes for C2.US panel
              dat <- fm$yHat # extracts predicted phenotypes from BGLR output
            dat <- dat[(length(C2.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.US panel
  # 9.11 Cor.test to compare predicted an observed values! 
          Pearson.Cor <- cor.test(dat, C2.US.Y$ISW)
          return(Pearson.Cor$estimate)# 0.70 average heritability of ISW
  # 9.12 Clean up R environment
      rm(list= c("C2.US.NA", "C2.US.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 9.13 Close loop and function
    }
    )
  # 9.14 unlist and convert to numeric so I can calculate average  
    C2.GRIN.Train.US.valid <- as.numeric(unlist(C2.GRIN.Train.US.valid))
  # 9.15  Add average to output dataframe and the standard deviation
    Output.data[1,2] <- mean(C2.GRIN.Train.US.valid)
    Output.data[1,3] <- sd(C2.GRIN.Train.US.valid)    
#####################################################################################      
# Section 10.0: C2.GRIN panel as training population for C2.SK panel  
  # 10.1 Start a loop to run 10 iterations
    C2.GRIN.Train.SK.valid <- lapply(1:10, function(x){
  # 10.2 Subset data frame to include C2.GRIN panel
    C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  # 10.3 Find those in the C2.SK panel
    C2.SK.Y <- filter(Pheno_Popped, Population == "C2.SK")
  # 10.4 Make object without C2.SK.Y ISW BLUEs
    C2.SK.NA <- C2.SK.Y
    C2.SK.NA$ISW <- NA
  # 10.5 Merge dataframes  
    Exp.Y <- rbind(C2.GRIN.Y, C2.SK.NA)
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
    dat <- dat[(length(C2.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.SK panel
  # 10.11 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, C2.SK.Y$ISW)
    return(Pearson.Cor$estimate)
  # 10.12 Clean up R environment
    rm(list= c("C2.SK.NA", "C2.SK.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))   
  # 10.13 Close loop and function
    }
    )
  # 10.14 unlist and convert to numeric so I can calculate average  
    C2.GRIN.Train.SK.valid <- as.numeric(unlist(C2.GRIN.Train.SK.valid))
  # 10.15  Add average to output dataframe and the standard deviation
    Output.data[2,2] <- mean(C2.GRIN.Train.SK.valid)
    Output.data[2,3] <- sd(C2.GRIN.Train.SK.valid)     
#####################################################################################      
# Section 11.0: C2.GRIN panel as training population for OH.GRIN panel 
  # 11.1 Start a loop to run 10 iterations
    C2.GRIN.Train.OHGRIN.valid <- lapply(1:10, function(x){  
  # 11.2 Subset data frame to include C2.GRIN panel
    C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")    
  # 11.3 Find those in the OH.GRIN panel
    OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")
  # 11.4 Make object without OH.GRIN.Y ISW BLUEs
    OH.GRIN.NA <- OH.GRIN.Y
    OH.GRIN.NA$ISW <- NA
  # 11.5 Merge dataframes  
    Exp.Y <- rbind(C2.GRIN.Y, OH.GRIN.NA)
  # 11.6 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 11.7 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 11.8 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 11.9 Run BGLR genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 11.10 Extract predicted phenotypes for OH.GRIN panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(C2.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.GRIN panel
  # 11.11 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, OH.GRIN.Y$ISW)
    return(Pearson.Cor$estimate)
    # 11.12 Clean up R environment
    rm(list= c("OH.GRIN.NA", "OH.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))       
  # 11.13 Close loop and function
    }
    )
  # 11.14 unlist and convert to numeric so I can calculate average  
    C2.GRIN.Train.OHGRIN.valid <- as.numeric(unlist(C2.GRIN.Train.OHGRIN.valid))
  # 11.15  Add average to output dataframe and the standard deviation
    Output.data[3,2] <- mean(C2.GRIN.Train.OHGRIN.valid)
    Output.data[3,3] <- sd(C2.GRIN.Train.OHGRIN.valid)   
#####################################################################################      
# Section 12.0: Leave-one-out Cross-validation within C2.GRIN panel 
  # 12.1 Subset data frame to include C2.GRIN panel
    C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")    
  # 12.2 Subset genotype data to just C2.GRIN
    Exp.G <- subset(Geno, taxa %in% C2.GRIN.Y$taxa)
  # 12.3 Remove "Taxa" column from genotype file
    Exp.G <- Exp.G[,-1]
  # 12.4 Subset phenotype data to just C2.GRIN
    Exp.Y <- C2.GRIN.Y
  # 12.5 Specify number of iterations to complete leave-one-out cross validation  
    LOOCV <- 1:length(Exp.Y$taxa)
  # 12.6 Creat list object to store estimated phenotypes  
    Loop.out.data <- list()
  # 12.7 start lapply loop for leave-one-out cross validation
    Loop.out.data <- lapply(LOOCV, function (i) { # LOOCV defined in 12.5
      ETA <- list(list(X=Exp.G, model='BL')) # specify parts of model
        Lp.Exp.Y <- Exp.Y # make object for phenotype data
          Lp.Exp.Y$ISW[i] <- NA # drop phenotype for LOOCV
          fm <- BGLR(y=Lp.Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
            # run model
          dat <- fm$yHat # get estimated phenotypes
        return(Loop.out.data <- dat[i]) # store estimated phenotype of left-out PI
      } 
    ) # repeat for number of individuals in populations
  # 12.8 Format LOOCV output into dataframe so cor.test can be done
    Hold.4.cor <- as.data.frame(unlist(Loop.out.data))
  # 12.9 Give informative column names
    colnames(Hold.4.cor)[1] <- "Out"
  # 12.10 Run Correlation  
    Pearson.Cor <- cor.test(Hold.4.cor$Out, Exp.Y$ISW)
  # 12.11 Calculate and store accuracy
    Output.data[4,2] <- Pearson.Cor$estimate
    Output.data[4,3] <- NA
  # Clean up R environment  
    rm(list= c("Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))      
#####################################################################################      
# Section 13.0 write output to csv file   
  # 13.1 give output more informative column names    
    colnames(Output.data) <- c("Populations", "Accuracy", "Stand.dev")
  # 13.2 write to csv file
    write.csv(x=Output.data, file="C2.GRIN.Trainer.csv", row.names=F)
  # 13.3 Prepare for next round analysis     
    rm(list= c("C2.SK.Y", "Output.data"))    
    Output.data <- matrix(nrow=4, ncol=3)
  # 13.4 Fill first column with Names of training and validation populations  
    Train.Valid.Names <- c("OH.GRIN Train and OH.US Valid",
                           "OH.GRIN Train and OH.SK Valid",
                           "OH.GRIN Train and C2.GRIN Valid",
                           "OH.GRIN LOOCV")
    Output.data[1:4,1] <- Train.Valid.Names    
#####################################################################################      
# Section 14.0: OH.GRIN panel as training population for OH.US panel
  # 14.1 Start a loop to run 10 iterations
    OH.GRIN.Train.OH.US.valid <- lapply(1:10, function(x){
  # 14.2 Filter the Phenotype file to find PIs assigned to OH.GRIN panel
    OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")
  # 14.3 Find those in the OH.US panel
    OH.US.Y <- filter(Pheno_Popped, Population == "OH.US")
  # 14.4 Make object without OH.US ISW BLUEs
    OH.US.NA <- OH.US.Y
    OH.US.NA$ISW <- NA
  # 14.5 Merge dataframes  
    Exp.Y <- rbind(OH.GRIN.Y, OH.US.NA)
  # 14.6 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 14.7 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 14.8 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 14.9 Run BGLR genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 14.10 Extract predicted phenotypes for OH.US panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(OH.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.US panel
  # 14.11 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, OH.US.Y$ISW)
    return(Pearson.Cor$estimate)
    # 14.12 Clean up R environment
    rm(list= c("OH.US.NA", "OH.US.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 14.3 Close loop and function  
    }
    )
  # 14.14 unlist and convert to numeric so I can calculate average  
    OH.GRIN.Train.OH.US.valid <- as.numeric(unlist(OH.GRIN.Train.OH.US.valid))
  # 14.15 Add average to output dataframe and the standard deviation
    Output.data[1,2] <- mean(OH.GRIN.Train.OH.US.valid)
    Output.data[1,3] <- sd(OH.GRIN.Train.OH.US.valid)  
#####################################################################################      
# Section 15.0: OH.GRIN panel as training population for OH.SK panel  
  # 15.1 State a loop for 10 iterations
    OH.GRIN.Train.OH.SK.valid <- lapply(1:10, function(x){
  # 15.2 Filter the Phenotype file to find PIs assigned to OH.GRIN panel
    OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")    
  # 15.3 Find those in the OH.SK panel
    OH.SK.Y <- filter(Pheno_Popped, Population == "OH.SK")
  # 15.4 Make object without OH.SK ISW BLUEs
    OH.SK.NA <- OH.SK.Y
    OH.SK.NA$ISW <- NA
  # 15.5 Merge dataframes  
    Exp.Y <- rbind(OH.GRIN.Y, OH.SK.NA)
  # 15.6 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 15.7 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 15.8 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 15.9 Run BGLR genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 15.10 Extract predicted phenotypes for OH.SK panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(OH.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the OH.SK panel
  # 15.11 Cor.test to compare predicted an observed values!   
    Pearson.Cor <- cor.test(dat, OH.SK.Y$ISW)
    return(Pearson.Cor$estimate)
    # 15.12 Clean up R environment
    rm(list= c("OH.SK.NA", "OH.SK.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 15.3 Close loop and function  
    }
    )
  # 15.14 unlist and convert to numeric so I can calculate average  
    OH.GRIN.Train.OH.SK.valid <- as.numeric(unlist(OH.GRIN.Train.OH.SK.valid))
  # 15.15 Add average to output dataframe and the standard deviation
    Output.data[2,2] <- mean(OH.GRIN.Train.OH.SK.valid)
    Output.data[2,3] <- sd(OH.GRIN.Train.OH.SK.valid)  
#####################################################################################      
# Section 16.0: OH.GRIN panel as training population for C2.GRIN panel  
  # 16.1 State loop for 10 iterations
    OH.GRIN.Train.C2SK.valid  <- lapply(1:10, function(x){
  # 16.2 Filter the Phenotype file to find PIs assigned to OH.GRIN panel
    OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")    
  # 16.3 Find those in the C2.GRIN panel
    C2.GRIN.Y <- filter(Pheno_Popped, Population == "C2.GRIN")
  # 16.4 Make object without C2.GRIN.Y ISW BLUEs
    C2.GRIN.NA <- C2.GRIN.Y
    C2.GRIN.NA$ISW <- NA
  # 16.5 Merge dataframes  
    Exp.Y <- rbind(OH.GRIN.Y, C2.GRIN.NA)
  # 16.6 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 16.7 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 16.8 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 16.9 Run BGLR genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 16.10 Extract predicted phenotypes for C2.GRIN panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(OH.GRIN.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.GRIN panel
  # 16.11 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, C2.GRIN.Y$ISW)
    return(Pearson.Cor$estimate)
    # 16.12 Clean up R environment
    rm(list= c("C2.GRIN.NA", "C2.GRIN.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
  # 16.13 Close loop and function
    }
    )
  # 16.14 unlist and convert to numeric so I can calculate average  
    OH.GRIN.Train.C2SK.valid <- as.numeric(unlist(OH.GRIN.Train.C2SK.valid))
  # 16.15 Add average to output dataframe and the standard deviation
    Output.data[3,2] <- mean(OH.GRIN.Train.C2SK.valid)
    Output.data[3,3] <- sd(OH.GRIN.Train.C2SK.valid)      
#####################################################################################      
# Section 17.0: Leave-one-out Cross-validation within OH.GRIN panel 
  # 17.1 State the validation and experimental lines for output
    OH.GRIN.Y <- filter(Pheno_Popped, Population == "OH.GRIN")    
  # 17.2 Subset genotype data to just OH.GRUN
    Exp.G <- subset(Geno, taxa %in% OH.GRIN.Y$taxa)
  # 17.3 Remove "Taxa" column from genotype file
    Exp.G <- Exp.G[,-1]
  # 17.4 Subset phenotype data to just OH.GRIN
    Exp.Y <- OH.GRIN.Y
  # 17.5 Specify number of iterations to complete leave-one-out cross validation  
    LOOCV <- 1:length(Exp.Y$taxa)
  # 17.6 Creat list object to store estimated phenotypes  
    Loop.out.data <- list()
  # 17.7 start lapply loop for leave-one-out cross validation
    Loop.out.data <- lapply(LOOCV, function (i) { # LOOCV defined in 12.5
      ETA <- list(list(X=Exp.G, model='BL')) # specify parts of model
        Lp.Exp.Y <- Exp.Y # make object for phenotype data
          Lp.Exp.Y$ISW[i] <- NA # drop phenotype for LOOCV
          fm <- BGLR(y=Lp.Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE) # run model
          dat <- fm$yHat # get estimated phenotypes
        return(Loop.out.data <- dat[i]) # store estimated phenotype of left-out PI
      } 
    ) # repeat for number of individuals in populations
  # 17.8 Format LOOCV output into dataframe so cor.test can be done
    Hold.4.cor <- as.data.frame(unlist(Loop.out.data))
  # 17.9 Give informative column names
    colnames(Hold.4.cor)[1] <- "Out"
  # 17.10 Run Correlation  
    Pearson.Cor <- cor.test(Hold.4.cor$Out, Exp.Y$ISW)
  # 17.11 Calculate and store accuracy
    Output.data[4,2] <- Pearson.Cor$estimate
    Output.data[4,3] <- NA
#####################################################################################      
# Section 18.0 write output to csv file   
  # 18.1 give output more informative column names    
    colnames(Output.data) <- c("Populations", "Accuracy", "Standard Deviation")
  # 18.2 write to csv file
    write.csv(x=Output.data, file="OH.GRIN.Trainer.csv", row.names=F)      
# :)    