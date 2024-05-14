#!/bin/sh

echo $(date)

input_dir=/data1/cwas_input/SMC_input
output_dir=/data1/cwas_output/SMC_output.annotation_v6.1

cwas effective_num_test \
-i $output_dir/correlation/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.correlation_matrix.all.zarr \
-c_count $output_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_counts.txt \
-o_dir $output_dir/eff_num_test \
-s $input_dir/pheno_SMC_batch1-7_RCL40_visual_1559.tsv \
-ef

echo $(date)
echo "Done"

