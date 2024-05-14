## Case-control comparison based on Amyloid beta

# 1. CWAS annotation, categorizaiton, burden test
Run_CWAS_01_amyloid_beta_annotation_burden.sh

# 2. CWAS feature selection
Run_CWAS_02_amyloid_beta_feature_selection.sh

## After feature selection, subset analysis should be conducted on categories featuring positive R2 annotations

# 3. CWAS risk score after feature selection
Run_CWAS_03_amyloid_beta_RS_after_feature_selection.sh

# 4. CWAS permutation test
Run_CWAS_04_amyloid_beta_permutation.sh

# 5. CWAS burden shift
Run_CWAS_05_amyloid_beta_burdenshift.sh

# 6. CWAS correlation
Run_CWAS_06_amyloid_beta_correlation.sh

# 7. CWAS effective number of test
Run_CWAS_07_amyloid_beta_eff_num_test.sh

# 8. CWAS DAWN
Run_CWAS_08_amyloid_beta_DAWN.step1.sh
Run_CWAS_08_amyloid_beta_DAWN.step2.sh

