# Results:Genetic distance and diversity of whole experiments and combined panels
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
Sample.Geno.Pop <- data.frame(Plant.Intros=row.names(Sample.Geno.DF), Sample.Geno.DF)
# 4.3 One PI is duplicated in the SK.2016 and OH.2018 experiment
  # Remove second-duplicated PI
   Genet.An.Geno[-394,]
  # Match column name between Exp.panel.IDs and Sample.Geno.pop    
colnames(Sample.Geno.Pop)[1] <- "Taxa"
# 4.4 Add population to dataframe so populations can be distringuished
Genet.An.Geno <- merge(Sample.Geno.Pop, Exp.panel.IDs, by="Taxa")

#####################################################################################
# Section 5.0 Reference where each population is. 
# 5.1 Find PIs that are C2.US
C2.US <- which(Genet.An.Geno[,1002]=="C2.US")
# 5.2 Find PIs that are C2.GRIN
C2.GRIN <- which(Genet.An.Geno[,1002]=="C2.GRIN")
# 5.3 Find PIs that are C2.SK
C2.SK <- which(Genet.An.Geno[,1002]=="C2.SK")
# 5.4 Find PIs that are OH.US
OH.US <- which(Genet.An.Geno[,1002]=="OH.US")
# 5.5 Find PIs that are OH.GRIN
OH.GRIN <- which(Genet.An.Geno[,1002]=="OH.GRIN")
# 5.6 Find PIs that are OH.GRIN
OH.SK <- which(Genet.An.Geno[,1002]=="OH.SK")

#####################################################################################
# Section 6.0 Combine US and GRIN populations for Fst & He
# 6.1 Create temporary object to store data. 
TmpGene.An.Geno <- Genet.An.Geno
# 6.2 Change GRIN & US to same populations
TmpGene.An.Geno[C2.US,1002] <- "Tmp1"
TmpGene.An.Geno[C2.GRIN,1002] <- "Tmp1"
TmpGene.An.Geno[OH.US,1002] <- "Tmp2"
TmpGene.An.Geno[OH.GRIN,1002] <- "Tmp2"
# 6.1 Identify PIs included
ind <- as.character(TmpGene.An.Geno$Taxa)
# 6.2 Identify Populations
population <- as.character(TmpGene.An.Geno$Population)
# 6.3 Identify the genotype scores
locus <- TmpGene.An.Geno[,2:1001]
# 6.4 Make objected for adegenet and heifstat
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
# 6.5 calculate Fst
Fst.Res <- genet.dist(Mydata1, method = "WC84")
# WC84 = Weir and Cockerham 1984 publication
# 6.5 Extract results frm Fst.Res object
Fst.Res
# 6.6 Subset data to only include C2.US and C2.GRIN
Temp.G <- subset(TmpGene.An.Geno, Population == "Tmp1")  
# 6.7 Redifine individuals
ind <- as.character(Temp.G$Taxa)
# 6.8 Redifine populations
population <- as.character(Temp.G$Population)
# 6.9 Redifine locus
locus <- Temp.G[,2:1001]
# 6.10 Calculate He 
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)

# 6.11 Repeat for OH.US and OH.GRIN
Temp.G <- subset(Genet.An.Geno, Population == "Tmp2")   
ind <- as.character(Temp.G$Taxa)
population <- as.character(Temp.G$Population)
locus <- Temp.G[,2:1001]
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)

#####################################################################################
# Section 7.0 Combine US and SK populations for Fst & He
# 7.1  Change SK & US to same populations 
TmpGene.An.Geno <- Genet.An.Geno
TmpGene.An.Geno[C2.US,1002] <- "Tmp1"
TmpGene.An.Geno[C2.SK,1002] <- "Tmp1"
TmpGene.An.Geno[OH.US,1002] <- "Tmp2"
TmpGene.An.Geno[OH.SK,1002] <- "Tmp2"
# 7.2 Identify PIs included
ind <- as.character(TmpGene.An.Geno$Taxa)
# 7.3 Identify Populations
population <- as.character(TmpGene.An.Geno$Population)
# 7.4 Identify the genotype scores
locus <- TmpGene.An.Geno[,2:1001]
# 7.5 Make objected for adegenet and heifstat
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
# 7.6 calculate Fst
Fst.Res <- genet.dist(Mydata1, method = "WC84")
# WC84 = Weir and Cockerham 1984 publication
# 7.7 Extract results frm Fst.Res object
Fst.Res
# 7.8 Subset data to only include C2.US and C2.SK
Temp.G <- subset(Genet.An.Geno, Population == c("C2.US", "C2.SK"))  
# 7.9 Redifine individuals
ind <- as.character(Temp.G$Taxa)
# 7.10 Redifine populations
population <- as.character(Temp.G$Population)
# 7.11 Redifine locus 
locus <- Temp.G[,2:1001]
# 7.12 Calculate He 
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)

# 7.13 Repeat for OH.US and OH.SK
Temp.G <- subset(Genet.An.Geno, Population == c("OH.US", "OH.SK"))  
ind <- as.character(Temp.G$Taxa)
population <- as.character(Temp.G$Population)
locus <- Temp.G[,2:1001]
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)

#####################################################################################
# Section 8.0 Combine GRIN and SK populations for Fst & He
# 8.1  Change GRIN & SK to same populations 
TmpGene.An.Geno <- Genet.An.Geno
TmpGene.An.Geno[C2.GRIN,1002] <- "Tmp1"
TmpGene.An.Geno[C2.SK,1002] <- "Tmp1"
TmpGene.An.Geno[OH.GRIN,1002] <- "Tmp2"
TmpGene.An.Geno[OH.SK,1002] <- "Tmp2"
# 8.2 Identify PIs included
ind <- as.character(TmpGene.An.Geno$Taxa)
# 8.3 Identify Populations
population <- as.character(TmpGene.An.Geno$Population)
# 8.4 Identify the genotype scores
locus <- TmpGene.An.Geno[,2:1001]
# 8.5 Make objected for adegenet and heifstat
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
# 8.6 calculate Fst
Fst.Res <- genet.dist(Mydata1, method = "WC84")
# WC84 = Weir and Cockerham 1984 publication
# 8.7 Extract results frm Fst.Res object
Fst.Res
# 8.8 Subset data to only include C2.GRIN and C2.SK
Temp.G <- subset(Genet.An.Geno, Population == c("C2.SK", "C2.GRIN"))  
# 8.9 Redifine Indivduals
ind <- as.character(Temp.G$Taxa)
# 8.10 Redifine Population
population <- as.character(Temp.G$Population)
# 8.11 Redifine Locus
locus <- Temp.G[,2:1001]
# 8.12 Calculate He
TempG.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(TempG.Exp)
# 8.13 Repeat for OH.GRIN and OH.SK
Temp.G <- subset(Genet.An.Geno, Population == c("OH.SK", "OH.GRIN"))  
Temp.G$Population = "Combination"
ind <- as.character(Temp.G$Taxa)
population <- as.character(Temp.G$Population)
locus <- Temp.G[,2:1001]
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)

#####################################################################################
# Section 9.0 Calculate Fst and He for C2.2018 and OH.2018 Experiments
# 9.1  Change so all populations are the same. 
TmpGene.An.Geno <- Genet.An.Geno
TmpGene.An.Geno[C2.GRIN,1002] <- "Tmp1"
TmpGene.An.Geno[C2.SK,1002] <- "Tmp1"
TmpGene.An.Geno[C2.US,1002] <- "Tmp1"
TmpGene.An.Geno[OH.GRIN,1002] <- "Tmp2"
TmpGene.An.Geno[OH.SK,1002] <- "Tmp2"
TmpGene.An.Geno[OH.US,1002] <- "Tmp2"
# 9.2 Identify PIs included
ind <- as.character(TmpGene.An.Geno$Taxa)
# 9.3 Identify Populations
population <- as.character(TmpGene.An.Geno$Population)
# 9.4 Identify the genotype scores
locus <- TmpGene.An.Geno[,2:1001]
# 9.5 Make objected for adegenet and heifstat
Mydata1 <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
# 9.6 calculate Fst
Fst.Res <- genet.dist(Mydata1, method = "WC84")
# WC84 = Weir and Cockerham 1984 publication
# 9.7 Extract results frm Fst.Res object
Fst.Res
# 9.8 Subset data to only include C2.2018 data
Temp.G <- subset(Genet.An.Geno, Population == c("C2.SK", "C2.GRIN", "C2.SK"))
# 9.9 redifine individual
ind <- as.character(Temp.G$Taxa)
# 9.10 redifine population
population <- as.character(Temp.G$Population)
# 9.11 redifine locus
locus <- Temp.G[,2:1001]
# 9.12 Calculate He
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)
# 9.13 Repeat for OH.2018 data
Temp.G <- subset(Genet.An.Geno, Population == c("OH.SK", "OH.GRIN", "OH.SK"))  
ind <- as.character(Temp.G$Taxa)
population <- as.character(Temp.G$Population)
locus <- Temp.G[,2:1001]
Temp.G.Exp <- df2genind(locus, ploidy = 2, ind.names = ind, pop = population, sep="")
basic.stats(Temp.G.Exp)


