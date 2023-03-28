# Results:Genomic prediction across experiments to generate or validate models
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
    source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
    source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
  # 4.2 Create and move to directory to store the following output
    dir.create("GAPIT.Output/")
    setwd("GAPIT.Output/")
  # 4.3 Run GAPIT numeric conversion
    myGAPIT <- GAPIT(G=Geno, output.numerical=TRUE, file.path = "GAPIT.Output/" )
    # writes numeric file to GAPIT.Output/ directory
  # 4.4 read output files into R environment 
    Geno <- read_delim(file = "GAPIT.Genotype.Numerical.txt", col_names = TRUE, delim = "\t")
    Geno <-  Geno[-which(duplicated(Geno$taxa)=="TRUE"),]
  # 4.5 Move back to previous directory
    setwd("../")
#####################################################################################
# Section 5.0: Clean up R environment.
  # 5.1 identify R objects no longer needed
    Extras <- ls()
    Extras <- Extras[c(1:121,124:134)]
  # 5.2 Remove those objects
    rm(list=Extras)
    rm(Extras)
#####################################################################################
# Section 6.0: Load phenotypic data to R environment 
  # 6.1 Read in phenotype data for all three experiments 
    C2.2018.Y <- read_delim(file="Input.phenotype/C2.2018.pheno.csv", col_names = T, delim = ",")
    OH.2018.Y <- read_delim(file="Input.phenotype/OH.2018.pheno.csv", col_names = T, delim = ",")    
    SK.2016.Y <- read_delim(file="Input.phenotype/SK.2016.pheno.csv", col_names = T, delim = ",")
  # 6.2 Combine all three dataframes
    Pheno <- rbind(C2.2018.Y,OH.2018.Y,SK.2016.Y)
  # 6.3 Match Phenotypic data to Genotypic data
    Geno <- Geno[order(Geno$taxa),]
    Pheno <- Pheno[order(Pheno$taxa),]
  # 6.4 Remove Extra R objects
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
    nIter=20000
  # 8.2 define burn in value
    burnIn=5000
  # 8.3 Make output matrix to store data
    Output.data <- matrix(nrow=2, ncol=2)
#####################################################################################      
# Section 9.0: C2.2018 Exp as training population for 2016.SK panel
  # 9.1 State the validation and experimental lines for output
    Output.data[1,1] <- "C2.2018 Train and 2016.SK Valid"
  # 9.2 Filter the Phenotype file to find PIs assigned to C2.2018 panel
    C2.2018.Y <- filter(Pheno_Popped, Population == "C2.SK" | Population ==  "C2.US"| Population == "C2.GRIN")
  # 9.3 Find those in the 2016.SK panel
    R.SK.Y <- filter(Pheno_Popped, Population == "R.SK")
  # 9.4 Make object without 2016.SK ISW BLUEs
    R.SK.NA <- R.SK.Y
    R.SK.NA$ISW <- NA
  # 9.5 Merge dataframes  
    Exp.Y <- rbind(C2.2018.Y, R.SK.NA)
  # 9.6 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 9.7 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 9.8 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 9.9 Run BGLR genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 9.10 Extract predicted phenotypes for 2016.SK panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(C2.2018.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.US panel
  # 9.11 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, R.SK.Y$ISW)
    Output.data[1,2] <- Pearson.Cor$estimate 
  # 9.12 Clean up R environment
    rm(list= c("C2.2018.Y", "Exp.G", "Exp.Y", "fm", "ETA", "dat", "Pearson.Cor"))
#####################################################################################      
# Section 10.0: OH.2018 Exp as training population for 2016.SK panel
  # 10.1 State the validation and experimental lines for output
    Output.data[2,1] <- "OH.2018 Train and 2016.SK Valid"
  # 10.2 Filter the Phenotype file to find PIs assigned to OH.2018 panel
    OH.2018.Y <- filter(Pheno_Popped, Population == "OH.SK" |Population == "OH.US" |Population == "OH.GRIN")
  # 10.3 Merge dataframes  
    Exp.Y <- rbind(OH.2018.Y, R.SK.NA)
  # 10.4 Remove Population Assignment
    Exp.Y <- Exp.Y[,-3]
  # 10.5 Make genotype file
    Exp.G <- subset(Geno, taxa %in% Exp.Y$taxa)
  # 10.6 Match order between dataframes!
    Exp.G <- Exp.G[match(Exp.Y$taxa, Exp.G$taxa),]    
  # 10.7 Run BGLR genomic Prediction
    ETA<-list(list(X=Exp.G[,-1], model='BL'))
    fm <- BGLR(y=Exp.Y$ISW, ETA=ETA, nIter=nIter, burnIn=burnIn, verbose=FALSE)
  # 10.10 Extract predicted phenotypes for 2016.SK panel
    dat <- fm$yHat # extracts predicted phenotypes from BGLR output
    dat <- dat[(length(OH.2018.Y$taxa)+1):length(dat)] # subset to only consider those for the C2.US panel
  # 10.11 Cor.test to compare predicted an observed values! 
    Pearson.Cor <- cor.test(dat, R.SK.Y$ISW)
    Output.data[2,2] <- Pearson.Cor$estimate 
#####################################################################################      
# Section 11.0 write output to csv file   
  # 11.1 give output more informative column names    
      colnames(Output.data) <- c("Populations", "Accuracy")
  # 11.2 write to csv file
      write.csv(x=Output.data, file="Cross.Exp.csv", row.names=F)      
# :)    
                    