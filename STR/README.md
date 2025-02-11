# AD Short tandem repeat analysis

We partially followed the workflow of [Guo et al., 2023] (GitHub: https://github.com/mhguo1/AD_STR) for STR analysis.

Specifically:
- 01.EH_preprocessing.R followed a similar input parsing and coverage extraction logic as parse_eh_coverage.R and parse_eh_genotypes.R, maintaining the same data processing pipeline while adjusting for dataset-specific formats.
- 03.EH_model1_model2.R was structured similarly to single_str_association_analysis.R and fisher_single_str_test.R, implementing analogous statistical tests but adapted to our dataset.
- 04.EH_DBSCAN_outlier.R implemented the same DBSCAN clustering approach as in run_dbscan.R and parse_dbscan.R, using identical distance metrics and parameter settings.

Additionally:
- 02.EH_QC.R was newly implemented for quality control based on coverage.
- 05.EH_model3.R was developed for STR analysis based on outlier counts.

For STR genotyping, we used ExpansionHunter v5 (https://github.com/Illumina/ExpansionHunter). The command script used is shown below, with the --variant-catalog option referencing a pre-established file from a previous study (https://github.com/mhguo1/AD_STR/tree/main/STR_genotyping).
- EH.sh

Below is a summary of our code.
1. Processing data from ExpansionHunter output vcf (01.EH_preprocessing.R)
2. Quality control using STR coverage (02.EH_QC.R)
3. STR analysis based on model 1, model 2 (03.EH_model1_model2.R)
4. STR outlier calling by DBSCAN for model 3 (04.EH_DBSCAN_outlier.R)
5. STR analysis based on model 3 (05.EH_model3.R)
