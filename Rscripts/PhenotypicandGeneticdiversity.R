# Results:Phenotypic and Genetic diversity between PI panels
# Author: William Rolling
# Corresponding Author: Leah McHale (mchale.21@osu.edu)
#####################################################################################
# Section 1.0: Prepare R environment
  # 1.1 Clear R environment
    rm(list=ls()) 
  # 1.2 List packages need to complete analysis
    Package.List <- c("adegenet",
                      "hierfstat",
                      "data.table",
                      "bigmemory",
                      "biganalytics",
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
  # 1.5 Some random number generation below; setting seed to make it repeatable
    set.seed(1122)  
#####################################################################################
# Section 2.0: Load genotypic data into R environment    
  # 2.1 Load genotype files for three experiments
    SK.2016 <- read_csv("Input.data/SK.2016.hmp.csv", col_names=TRUE)
    C2.2018 <- read_csv("Input.data/C2.2018.hmp.csv", col_names=TRUE)
    colnames(C2.2018)[1] <- "rs"
    OH.2018 <- read_csv("Input.data/OH.2018.hmp.csv", col_names=TRUE)
    colnames(OH.2018)[1] <- "rs"
  # 2.3 Match SNPs across all three populations in SK.2016 Experiment
    SK.2016.formatted <- filter(SK.2016, rs %in% C2.2018$rs) 
    SK.2016.formatted <- filter(SK.2016.formatted, rs %in% OH.2018$rs) 
  # 2.4 Remove hmp file format extra columns
    SK.2016.formatted <- select(SK.2016.formatted, -(2:11))
  # 2.5 order column in ascending order on SNP name
    SK.2016.formatted <- SK.2016.formatted[order(SK.2016.formatted$rs),]
  # 2.5 Repeat steps 2.3 and 2.4 for OH.2018 and C2.2018 experiments  
    C2.2018.formatted <- filter(C2.2018,rs %in% SK.2016.formatted$rs) 
    C2.2018.formatted <- C2.2018.formatted[order(C2.2018.formatted$rs),]
    OH.2018.formatted <- filter(OH.2018,rs %in% SK.2016.formatted$rs) 
    OH.2018.formatted <- OH.2018.formatted[order(OH.2018.formatted$rs),]
     # Can remove first 12 columns because SNP names match SK dataframe!
    C2.2018.formatted <- select(C2.2018.formatted, -(1:11))
    OH.2018.formatted <- select(OH.2018.formatted, -(1:11))
#####################################################################################
# Section 3.0: Combine all genotype files and format     
  # 3.1 bind all genotype data vertically
    Geno <- cbind(SK.2016.formatted,C2.2018.formatted,OH.2018.formatted)
  # 3.2 Sample 1000 markers
        # Only a subset of markers are needed for this analysis. 
        # Sample 1000 markers here to speed up the following
    Sub.sample <- sample(2:length(Geno$rs), size=1000, replace=FALSE) %>%
                    sort()
    Sample.Geno <- Geno[Sub.sample,]
  # 3.3 Format genotype object by transposing as a dataframe   
    Sample.Geno.DF <- as.data.frame(Sample.Geno) %>%
                            t()
    colnames(Sample.Geno.DF) <- Sample.Geno.DF[1,]
    Sample.Geno.DF <- Sample.Geno.DF[-1,] 
#####################################################################################
# Section 4.0: Load Population assignments     
  # 4.1 Load  population assignment file created during Rolling et al. 2020a
      Exp.panel.IDs <- read_csv("Input.data/Experiment_Panel.IDs.csv", col_names = TRUE)
  # 4.2 Make sure Taxa matches across genotype and population IDs
      Sample.Geno.Pop <- data.frame(Plant.Intros=row.names(Sample.Geno.DF),Sample.Geno.DF)
  # 4.3 One PI is duplicated in the SK.2016 and OH.2018 experiment
      # Remove second-duplicated PI
      Removes <- which(duplicated(Sample.Geno.Pop$Plant.Intros) ==TRUE)
      Sample.Geno.Pop <- Sample.Geno.Pop[-Removes,]
  # Match column name between Exp.panel.IDs and Sample.Geno.pop    
      colnames(Sample.Geno.Pop)[1] <- "Taxa"
  # 4.4 Add population to dataframe so populations can be distringuished
      Genet.An.Geno <- merge(Sample.Geno.Pop, Exp.panel.IDs, by="Taxa")
#####################################################################################
# Section 5.0: Define Objects for Fst Analysis:
  # 5.1 Identify PIs included
      ind <- as.character(Genet.An.Geno$Taxa)
  # 5.2 Identify Populations
      population <- as.character(Genet.An.Geno$Population)
  # 5.3 Identify the genotype scores
      locus <- Genet.An.Geno[,2:1001]
  # 5.4 Make objected for adegenet and heifstat
      Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
#####################################################################################
# Section 6.0: Calculate Fst and make dpac plto
  # 6.1 calculate Fst
      Fst.Res <- genet.dist(Mydata1, method = "WC84")
      # WC84 = Weir and Cockerham 1984 publication
  # 6.2 Extract results frm Fst.Res object
      Fst.Res
        # I copied results from R console to excel
  # 6.3 create scatter plot values
      dapc1 <- dapc(Mydata1, population, n.pca=20, n.da=2) 
        # Have population assignements in population column
        # I tried n.pca = 20, 40, 100 results are similar
        # I tried n.da = 2, 6, 50 results are similar
#####################################################################################      
# 7.0 Make scatter plot similar to Figure 1a
  # 7.1 Define colors for different panels of PIs    
      myCol <- c("Grey","Grey","Red","Red","Black","Black", "darkgreen")
  # 7.2 Define shapes for different panels of PIs    
      myPch <- c(15:20)
  # 7.3 Prepare to write scatter plot to tiff figure 
      tiff(filename = "Figure1a.tiff", bg="white")
  # 7.4 Draw scatter plot      
      scatter(dapc1, posi.da="NULL", bg="white",
              pch=myPch, cstar=0, col=myCol, scree.pca=TRUE,
              posi.pca="bottomleft",label = NULL)
      # Throws an error but it still makes the plot I want. 
  # 7.5 shut down "device" to write to tiff file
      dev.off()
#####################################################################################      
# 8.0 Calculate within population diversity... commented example with SK.2016 experiment
  # 8.1 subset based on Population ID in Genet.An.Geno
      SK.2016.G <- subset(Genet.An.Geno, Population == "R.SK")  
  # 8.2 Identify PIs included
      ind <- as.character(SK.2016.G$Taxa)
  # 8.3 Identify Populations
      population <- as.character(SK.2016.G$Population)
  # 8.4 Identify the genotype scores
      locus <- SK.2016.G[,2:1001]
  # 8.5 Make objected for adegenet and heifstat
      SK.2016.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
  # 8.6 Extract Hs from "basic.stats" 
      basic.stats(SK.2016.Exp)
  # Hs = 0.25, Ho = 0.004, Ht = 0.2469
      
#####################################################################################           
  # 9.0 repeat for other experiments and panels
    # 9.1 C2.SK
      C2.SK.G <- subset(Genet.An.Geno, Population == "C2.SK")  
      ind <- as.character(C2.SK.G$Taxa)
      population <- as.character(C2.SK.G$Population)
      locus <- C2.SK.G[,2:1001]
      C2.SK.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
      basic.stats(C2.SK.Exp)
      # Hs = 0.23, Ho = 0.004, Ht = 0.23
    # 9.2 C2.US
      C2.US.G <- subset(Genet.An.Geno, Population == "C2.US")  
      ind <- as.character(C2.US.G$Taxa)
      population <- as.character(C2.US.G$Population)
      locus <- C2.US.G[,2:1001]
      C2.US.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
      basic.stats(C2.US.Exp)
      # Hs = 0.19, Ho = 0.008, Ht = 0.19  
    # 9.3 C2.US
      C2.GRIN.G <- subset(Genet.An.Geno, Population == "C2.GRIN")  
      ind <- as.character(C2.GRIN.G$Taxa)
      population <- as.character(C2.GRIN.G$Population)
      locus <- C2.GRIN.G[,2:1001]
      C2.GRIN.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
      basic.stats(C2.GRIN.Exp)
      # Hs = 0.39, Ho = 0.007, Ht = 0.39    
    # 9.4 OH.SK
      OH.SK.G <- subset(Genet.An.Geno, Population == "OH.SK")  
      ind <- as.character(OH.SK.G$Taxa)
      population <- as.character(OH.SK.G$Population)
      locus <- OH.SK.G[,2:1001]
      OH.SK.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
      basic.stats(OH.SK.Exp)
      # Hs = 0.23, Ho = 0.003, Ht = 0.23     
    # 9.5 OH.US
      OH.US.G <- subset(Genet.An.Geno, Population == "OH.US")  
      ind <- as.character(OH.US.G$Taxa)
      population <- as.character(OH.US.G$Population)
      locus <- OH.US.G[,2:1001]
      OH.US.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
      basic.stats(OH.US.Exp)
      # Hs = 0.19, Ho = 0.007, Ht = 0.39
    # 9.6 OH.GRIN
      OH.US.G <- subset(Genet.An.Geno, Population == "OH.US")  
      ind <- as.character(OH.US.G$Taxa)
      population <- as.character(OH.US.G$Population)
      locus <- OH.US.G[,2:1001]
      OH.US.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
      basic.stats(OH.US.Exp)  
    # Hs = 0.39, Ho = 0.007, Ht = 0.39    
      
