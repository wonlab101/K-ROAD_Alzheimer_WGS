#!/bin/sh

echo $(date)

input_dir=/input
output_dir=/output

# skip annotation (same as before)
cwas categorization -i $output_dir/annotation/CWAS_input.annotated.vcf.gz \
-o_dir $output_dir/categorization/ \
-p 40

cwas binomial_test \
-i $output_dir/categorization/CWAS_input.categorization_result.zarr \
-o_dir $output_dir/burden_test/ \
-s $input_dir/phenotype_table.tsv \
--use_n_carrier


echo $(date)
echo "Done"
