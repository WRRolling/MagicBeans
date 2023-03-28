# Results:Phenotypic and Genetic diversity between PI panels
# Author: William Rolling
# Corresponding Author: Leah McHale (mchale.21@osu.edu)
#####################################################################################
# Section 1.0: Prepare R environment
  # 1.1 Clear R environment
    rm(list=ls()) 
  # 1.2 Setworking directory
  # example for someone who downloads data to folder in desktop <---- Set pathway to folder here. 
    setwd("~/../Desktop/Supplemental Data 1/")
#####################################################################################
# Section 2.0: Load data into R environment    
  # 2.0 Load data
    Cor.dat <- read.csv("Input.data/Correlation.Fst.Accuracy.csv", head=T)
  # 2.1 Run correlation
    Res <- cor.test(Cor.dat$Fst, Cor.dat$Accuracy)
  # 2.2 Look at results 
    Res
