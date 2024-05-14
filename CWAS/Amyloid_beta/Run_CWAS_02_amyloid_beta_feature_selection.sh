#!/bin/sh

echo $(date)

input_dir=/data2/cwas_paper/cwas_input/SMC_input
output_dir=/data2/cwas_paper/cwas_output/SMC_output.annotation_v6.1

cwas risk_score -i $output_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.categorization_result.zarr/ \
-s $input_dir/pheno_SMC_batch1-7_RCL40_visual_1559.tsv \
-o_dir $output_dir/feature_selection/ \
-c_info $output_dir/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.category_info.txt \
--domain_list noncoding \
--use_n_carrier \
-thr 2 \
-n_reg 10 \
-f 5 \
-p 40 \
-tf 0.7 \
--predict_only \
-t 20240416.seed42.tf70 \
--do_each_one \
-fs_group gene_set,functional_annotation,functional_score

echo $(date)
echo "Done"