# Results:Genetic distance and diversity of whole experiments and combined panels
# Author: William Rolling
# Corresponding Author: Leah McHale (mchale.21@osu.edu)
#####################################################################################
# Section 1.0: Prepare R environment
# 1.1 Clear R environment
rm(list=ls()) 
# 1.2 List packages need to complete analysis
Package.List <- c("lme4",
                  "lmerTest",
                  "tidyverse",
                  "coefplot",
                  "ggpubr",
                  "cowplot",
                  "sjPlot",
                  "sjmisc",
                  "effects",
                  "sjstats",
                  "MuMIn")

# 1.3 Install and/or load R packages
# Script from: https://vbaliga.github.io
package.check <- lapply(
  Package.List ,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})
# 1.4 Load Data
dat <- read_csv("../Mod.dat.csv", col_names = T)
colnames(dat)[8] <- "Accuracy"
#############################################################
# Section 2 Load Phenotypic Data
# 2.1 Load phenotypes from three separate files
Pheno.2016 <- read_csv("../Input.phenotype/SK.2016.pheno.csv", col_names = T)
Pheno.OH <- read_csv("../Input.phenotype/OH.2018.pheno.csv", col_names = T)
Pheno.C2 <- read_csv("../Input.phenotype/C2.2018.pheno.csv", col_names = T)
# 2.2 Load Panel IDs
My.Panels <- read_csv("../Input.data/Experiment_Panel.IDs.csv")
# 2.3 Create Object to store output
Phenotype.Var <- matrix(nrow=9, ncol=2)
# 2.4 Combine Phenotypic data with the Model.data
Pheno.C2.Panel <- merge(Pheno.C2, My.Panels, 
                        by.x="taxa", by.y="Taxa")
# 2.5 Calculate Standard Deviation in C2.SK panel
  # Find individuals in C2. SK panel
  ref.no <- which(Pheno.C2.Panel$Population =="C2.SK")
  # Subset data to only include those in C2 Panel
  St.calc <- Pheno.C2.Panel[ref.no,]
  # Name row for output data
  Phenotype.Var[1,1] <- "C2.SK"  
  # Calculate standard deviation for C2.SK panel
  Phenotype.Var[1,2] <- sd(St.calc$ISW, na.rm = T) 
# 2.6 Repeat Standard Deviation Calculation for C2.GRIN panel 
  ref.no <- which(Pheno.C2.Panel$Population =="C2.GRIN")
  St.calc <- Pheno.C2.Panel[ref.no,]
  Phenotype.Var[2,1] <- "C2.GRIN"  
  Phenotype.Var[2,2] <- sd(St.calc$ISW, na.rm = T)
# 2.7 Repeat Standard Deviation Calculation for C2.US panel
  ref.no <- which(Pheno.C2.Panel$Population =="C2.US")
  St.calc <- Pheno.C2.Panel[ref.no,]
  Phenotype.Var[3,1] <- "C2.US"  
  Phenotype.Var[3,2] <- sd(St.calc$ISW, na.rm = T) 
# 2.8 Standard Deviation Calculation for all PIs in the 2018 C2 Population
  Phenotype.Var[4,1] <- "2018-C2"  
  Phenotype.Var[4,2] <- sd(Pheno.C2.Panel$ISW, na.rm = T) 
# 2.9 Repeat Standard Deviation Calculations for all OH Panels
  # OH.SK panel
  Pheno.OH.Panel <- merge(Pheno.OH, My.Panels, 
                          by.x="taxa", by.y="Taxa")  
  ref.no <- which(Pheno.OH.Panel$Population =="OH.SK")
  St.calc <- Pheno.OH.Panel[ref.no,]
  Phenotype.Var[5,1] <- "OH.SK"  
  Phenotype.Var[5,2] <- sd(St.calc$ISW, na.rm = T)   
  # OH.GRIN panel
  ref.no <- which(Pheno.OH.Panel$Population =="OH.GRIN")
  St.calc <- Pheno.OH.Panel[ref.no,]
  Phenotype.Var[6,1] <- "OH.GRIN"  
  Phenotype.Var[6,2] <- sd(St.calc$ISW, na.rm = T)   
  # OH.US panel
  ref.no <- which(Pheno.OH.Panel$Population =="OH.US")
  St.calc <- Pheno.OH.Panel[ref.no,]
  Phenotype.Var[7,1] <- "OH.US"  
  Phenotype.Var[7,2] <- sd(St.calc$ISW, na.rm = T) 
  # OH All
  Phenotype.Var[8,1] <- "2018-OH"  
  Phenotype.Var[8,2] <- sd(Pheno.OH.Panel$ISW, na.rm = T)  
# 2.10 Calculate Standard Deviation for the 2016-SK population
  Phenotype.Var[9,1] <- "2016-SK"
  Phenotype.Var[9,2] <- sd(Pheno.2016$ISW, na.rm=T)
# 2.11 Add Phenotypic Variation to Model.data
  Phenotype.Var <- as.data.frame(Phenotype.Var)
  Phenotype.Var$V2 <- as.numeric(Phenotype.Var$V2)
  dat.1 <- merge(dat, Phenotype.Var, 
                          by.x="Valid", by.y="V1")
#############################################################
# Section 3 Create Scatter Plots to consider the data!
# 3.1 Plot the Expected Heterozygosity of the Training Populations
HeTrain <- ggplot(dat.1, aes(x=`Train He`, y=Accuracy)) +
            geom_point() +
            xlab(label = "He Training Panel") +
            theme_classic()
# 3.2 Plot the Expected Heterozygosity of the Validation Population  
HeValid <- ggplot(dat.1, aes(x=`Valid He`, y=Accuracy)) +
  geom_point() +
  xlab(label = "He Validation Panel") +
  theme_classic()
# 3.3 Plot the Fst value between Training & Validation Population
Fst <- ggplot(dat.1, aes(x=Fst, y=Accuracy)) +
  geom_point() +
  xlab(label = "Fst Training & Validation Panels") +
  theme_classic()
# 3.4 Plot the number of individuals in the training population
NumTrain <- ggplot(dat.1, aes(`Num Train`, y=Accuracy)) +
  geom_point() +
  xlab(label = "Number of Individuals in Training Panel") +
  theme_classic()
# 3.5 Plot the number of individuals in the training population
NumValid <- ggplot(dat.1, aes(`Num Valid`, y=Accuracy)) +
  geom_point() +
  xlab(label = "Number of Individuals in Validation Panel") +
  theme_classic()
# 3.6 Plot the standard deviation in the validation population
SDValid <- ggplot(dat.1,aes(V2, y=Accuracy)) +
  geom_point() +
  xlab(label = "Standard Deviation of Phenotype in Validation Panel") +
  theme_classic()

png(filename = "Plotted.Data.png", bg="transparent", width=1000, height=500)
ggarrange(HeTrain, HeValid, Fst, NumTrain, 
                        NumValid, SDValid,  ncol=3, nrow=2)
dev.off()

#############################################################
# Section 4: Test which components are significantly related to Accuracy
# 4.1 Model w/ Fst
  # start w/ Fst because this (should) have largest effect on accuracy
  Mod.1 <- lmer(Accuracy ~ (1|Fst), data = dat.1) 
  summary(Mod.1) # Look at summary
  ranova(Mod.1) # Is Fst significant? 
# 4.2 Model w/ Fst & Validation Panel Genetic Diversity.
  # Next Largest effect?
  Mod.2 <- lmer(Accuracy ~ (1|Fst) + (1|`Valid He`), data = dat.1) 
  summary(Mod.2)
  ranova(Mod.2)
  anova(Mod.1, Mod.2) # see if models are significant different
  # Improvement to include Validation Panel - Marginal...
    # ValidHe significant component of the equation.
# 4.3 Model w/ Fst & Training Panel He. 
  Mod.3 <- lmer(Accuracy ~ (1|Fst) + (1|`Train He`), data = dat.1) 
  summary(Mod.3)
  ranova(Mod.3)
  anova(Mod.2, Mod.3)
    # Training He is not a significant component do not include!
# 4.4   Model with Fst Valid He and Phenotypic Standard Deviation
  Mod.4 <- lmer(Accuracy ~ (1|Fst) + (1|V2), data = dat.1)
  summary(Mod.4)
  ranova(Mod.4)
  anova(Mod.2, Mod.4)
  # Phenotypic Variation is significant
  # Use Mod.2 from now on!
# 4.5 Test Genetic Diversity & Phenotypic Diversity
  Mod.5 <- lmer(Accuracy ~1 + (1|Fst) + (1|V2) + (1|`Valid He`), data = dat.1)
  summary(Mod.5)
  ranova(Mod.5)
  anova(Mod.2, Mod.5)
  # Valid He is circa significance
  anova(Mod.2, Mod.5)
  # Use Mod.5 from now on!
  cor.test(dat.1$`Num Valid`, dat.1$`Valid He`)
  # negative correlation... Large panel is less diverse (SK-2016)


  
#############################################################
# Section 5: Test adding Fixed effects of number in training & validation panels.
  # 5.1 Test adding the number in the training panel
    Mod.5.1 <- lmer(Accuracy ~ 1 + `Num Train` + (1|Fst) + (1|V2) + (1|`Valid He`), data = dat.1) # An improvement
    summary(Mod.5.1)
    ranova(Mod.5.1)
    anova(Mod.5, Mod.5.1)
  # 5.2 Test adding the number in the validation panel
    Mod.5.2 <- lmer(Accuracy ~ 1 + `Num Valid` + (1|Fst) + (1|V2) + (1|`Valid He`), data = dat.1)
    anova(Mod.5.1, Mod.5.2)
  # 5.3 Test adding both the number of individuals in training and validation panels
    Mod.5.3 <- lmer(Accuracy ~ 1 + `Num Train`+ `Num Valid` + (1|Fst) + (1|`Valid He`), data = dat.1)  
    anova(Mod.5.3, Mod.5.1)
  # 5.4 Specify Best Model as Mod.final
       Mod.final <- Mod.5.1

#############################################################
# Section 6: Start Plotting & looking at variance explained
  # There are several iterations of different Fst & each panel. 
  # There is not a huge number of samples but can still check ICC
  # 5.5 measure ICC     
    performance::icc(Mod.5.1)
  
# Intraclass Correlation Coefficient
  #Adjusted ICC: 0.837
  #Unadjusted ICC: 0.760

# Measured R2
  r.squaredGLMM(Mod.final) #R2=0.09 Fixed & 0.85 Random
  r.squaredGLMM(Mod.5)    #R2=0.84
  r.squaredGLMM(Mod.2)  #R2=0.77
  r.squaredGLMM(Mod.1) #R2=0.60


dat.2 <- dat.1 
dat.2$model1_predicted <- predict(Mod.1)
dat.2$modelfinal_predicted <- predict(Mod.final)

# Plot Relationship between Fst & Accuracy 
png("linearMod.png",bg="transparent")
ggplot(data=dat.2,
       aes(x=Fst,
           y=Accuracy,))+
  geom_point(aes(colour=`Valid He`),
             size=(dat.2$`Num Train`/200),
             alpha=.7,
             position = "jitter") +
  scale_fill_gradient2(
    low="black",
    high="grey",
    midpoint = 0.2
  ) +
  geom_smooth(aes(y=model1_predicted), method=lm) +
  geom_smooth(aes(y=modelfinal_predicted), method=lm)+
  theme_classic() 
dev.off()




