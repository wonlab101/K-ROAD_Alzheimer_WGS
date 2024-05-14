#!/bin/sh

echo $(date)

input_dir=/data1/cwas_input/SMC_input
output_dir=/data1/cwas_output/SMC_output.annotation_v6.1

cwas correlation \
-i $output_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.categorization_result.zarr \
-cm sample \
-o_dir $output_dir/correlation \
-p 40 \
-c_info $output_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_info.txt

echo $(date)
echo "Done"

