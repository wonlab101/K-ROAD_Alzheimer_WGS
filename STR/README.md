## Short tandem repeat analysis
# Mainly followed the workflow of Guo et al., 2023

# Processing data from Expansionhunter output vcf
01.EH_preprocessing.R	

# Quality control
02.EH_QC.R		

# STR analysis based on model 1, model 2
03.EH_model1_model2.R	

# STR outlier calling by DBSCAN for model 3
04.EH_DBSCAN_outlier.R	

# STR analysis based on model 3
05.EH_model3.R	