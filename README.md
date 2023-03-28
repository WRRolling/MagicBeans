# MagicBeans

The files in these directories are the data used in the analyses associated with the manuscript " The effect of genetic distance and genetic diversity on genomic prediction accuracy for soybean quantitative disease resistance to Phytophthora sojae"
---------------------------------------------------------------------------------------
The data includes:
---------------------------------------------------------------------------------------
1) C2.2018.hmp.csv = Hapmap formatted genotypic data for germplasm panels C2.SK, C2.US, and C2.GRIN. 
2) OH.2018.hmp.csv = Hapmap formatted genotypic data for germplasm panels OH.SK, OH.US, and OH.GRIN.
3) SK.2016.hmp.csv = Hapmap formatted genotypic data for the 2016.SK germplasm panel.
4) C2.2018.pheno.csv = Inoculated shoot weight BLUE value for germplasm panels C2.SK, C2.US, and C2.GRIN. 
5) OH.2018.pheno.csv = Inoculated shoot weight BLUE value for ermplasm panels OH.SK, OH.US, and OH.GRIN. 
6) SK.2016.pheno.csv = Inoculated shoot weight BLUE value for the 2016.SK germplasm panel. 
7) Experiment_Panel.IDs.csv = Identification of which panels each plant introduction belongs to. 
7) Correlation.Fst.Accuracy.csv = Data used to measure the correlation between genetic distance and genomic prediction accuracy. 
8) Mod.dat.csv = Data used in regression analyses to identify which aspects of the training and validation panel have a significant effect on genomic prediction. 

The R scripts include:
-----------------------------------------------------------------------------------------
1) GeneticDiversityCombinedPanels.R - Analyses used to create Supplemental Figure 1, calculate Fst, and He
2) MakePhenotypeViolinPlots.R - Analyses used to create Supplemental Figure 1
3) SKtrainingpop.R - Analyses used to create Supplemental Table 1
4) UStrainingpop.R - Analyses used to create Supplemental Table 1
5) GRIN_germplasm_as_training_pop.R - Analyses used to create Supplemental Table 1
6) CombineTwopanelstoTrain.R - Analyses used to create Supplemental Table 1
7) USGRIN.Train2016SK.R - Analyses used to create Supplemental Table 1
8) 2016SKvalidate.R - Analyses used to create Supplemental Table 1
9) Correlation_Fst_Accuracy.R - Used to calculate correlation between Fst and genomic prediction accuracy. 
10) Regression.R - Used to complete regression analyses to identify which aspects of the training and validation panel have a significant effect on genomic prediction and create Figure 1, and provide data for supplemental table 2.  
11) AddRelated.OHSK2OHUS.R - Analyses used to create Table 2
12) AddRelated.OHUS2OHSk.R - Analyses used to create Table 2
13) AddRelated.US2GRIN.R - Analyses used to create Table 2
14) GeneticDiversityCombinedPanels.R - Analyses used to calculate He and Fst when two panels were combined to create training panel. Results in Supplemental Table 2. 
