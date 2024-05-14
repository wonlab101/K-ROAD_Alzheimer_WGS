#!/bin/sh

echo $(date)

input_dir=/data1/cwas_input/SMC_input
output_dir=/data1/cwas_output/SMC_output_969.annotation_v6.1

# skip annotation (same as before)
cwas categorization -i $output_dir/annotation/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.sorted.annotated.vcf.gz \
-o_dir $output_dir/categorization/ \
-p 40

cwas binomial_test \
-i $output_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.sorted.categorization_result.zarr \
-o_dir $output_dir/burden_test/ \
-s $input_dir/pheno_SMC_batch1-7_DXnoMCI_969_240422.tsv \
--use_n_carrier


echo $(date)
echo "Done"
