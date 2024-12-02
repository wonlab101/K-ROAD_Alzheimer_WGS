#!/bin/sh

echo $(date)

input_dir=/data2/cwas_paper/cwas_input/SMC_input
output_dir=/data2/cwas_paper/cwas_output/SMC_output.annotation_v6.1

cwas permutation_test -i $output_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.categorization_result.zarr \
-o_dir $output_dir/permutation_test/ \
-s $input_dir/pheno_SMC_batch1-7_RCL40_visual_1559.tsv \
-p 70 \
-n 10000 \
-b --use_n_carrier

echo $(date)
echo "Done"
