#!/bin/sh

echo $(date)

input_dir=/data2/cwas_paper/cwas_input/SMC_input
output_dir=/data2/cwas_paper/cwas_output/SMC_output.annotation_v6.1

cwas burden_shift \
-i $output_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.burden_test.txt \
-b $output_dir/permutation_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.binom_pvals.txt.gz \
-c_info $output_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_info.txt \
-c_count $output_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_counts.txt \
-o_dir $output_dir/burden_shift/ \
-c_cutoff 10

echo $(date)
echo "Done"
