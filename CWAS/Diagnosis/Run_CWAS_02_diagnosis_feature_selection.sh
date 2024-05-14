#!/bin/sh

echo $(date)

input_dir=/data1/cwas_input/SMC_input
output_dir=/data1/cwas_output/SMC_output_969.annotation_v6.1

cwas risk_score -i $output_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.sorted.categorization_result.zarr \
-s $input_dir/pheno_SMC_batch1-7_DXnoMCI_969_240422.tsv \
-o_dir $output_dir/feature_selection/ \
-c_info $output_dir/burden_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.sorted.category_info.txt \
--domain_list noncoding \
--use_n_carrier \
-t seed42.tf80 \
-thr 2 \
-tf 0.8 \
-p 40 \
--predict_only \
--do_each_one \
-fs_group gene_set,functional_annotation,functional_score

echo $(date)
echo "Done"