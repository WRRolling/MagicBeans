# Results:Distributions (Figure 1b)
# Author: William Rolling
# Corresponding Author: Leah McHale (mchale.21@osu.edu)
#####################################################################################
# Section 1.0: Prepare R environment
  # 1.1 Clear R environment
    rm(list=ls()) 
  # 1.2 List packages need to complete analysis
    Package.List <- c("ggpubr",
                  "hierfstat",
                  "cowplot",
                  "ggplot2",
                  "tidyverse")
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
  # 1.5 Set a color scheme for plot
    myColors <- c("red","white", "grey")

#####################################################################################
# Section 2.0: Make violin plot for C2.2018 data
  # 2.1 Load C2.2018 phenotypic data
    C2.dat <- read_csv("Input.phenotype/C2.2018.pheno.csv", col_names=T)
  # 2.2 Give columns informative names in C2.2018 dataset  
    colnames(C2.dat) <- c("Taxa", "ISW")
  # 2.3 Load  population assignment file created during Rolling et al. 2020a
    Exp.panel.IDs <- read_csv("Input.data/Experiment_Panel.IDs.csv", col_names = TRUE)
  # 2.4 Make sure Taxa matches across genotype and population IDs
    C2.dat.Pop <- merge(C2.dat, Exp.panel.IDs, by="Taxa")
  # 2.5 Prepare to write file plot to tiff file w/ white background
    tiff(filename = "Figure1bC2.tiff", bg="white")
  #2.6 ggplot to make violin plot 
    ggplot(C2.dat.Pop, aes(x=Population, y=ISW)) + # define data
      geom_violin(aes(fill=Population)) + # define type of plot as violin
      scale_fill_manual(values=myColors) + # add you color scheme from step 1.5
      theme_classic() + # I that theme_classic() keeps things uncluttered
      geom_boxplot(width=0.1) +  # Add a boxplot 
      ggtitle("2018-C2") + # Abbreviated name for plot: ISW=Inoculated Shoot Weight
      theme(legend.position = "none", # Aestictics: remove legend
        axis.title.x = element_blank(), # Aestictics; remove x axis label
        axis.title.y = element_blank(), # Aestictics: remove y axis label
        plot.title = element_text(hjust = 0.5)) # Aestictics: location of title
  # 2.7 turn of plotting to tiff file Figure 1b.   
    dev.off()
#####################################################################################
# Repeat for OH.2018 data 
    OH.dat <- read_csv("Input.phenotype/OH.2018.pheno.csv", col_names=T)
    colnames(OH.dat) <- c("Taxa", "ISW")
    OH.dat.Pop <- merge(OH.dat, Exp.panel.IDs, by="Taxa")
    # One extra step. There is one PI in common w/ 2016-SK dataset
    # Make sure this is defined as OH.SK here
    OH.dat.Pop$Population[which(OH.dat.Pop$Population == "R.SK")] <- "OH.SK"
    tiff(filename = "Figure1bOH.tiff", bg="white")
    ggplot(OH.dat.Pop, aes(x=Population, y=ISW)) +
      geom_violin(aes(fill=Population)) +
      scale_fill_manual(values=myColors) +
      theme_classic() +
      geom_boxplot(width=0.1) + 
      ggtitle("CRW") +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5)) 
    dev.off()
