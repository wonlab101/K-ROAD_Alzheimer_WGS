#!/bin/sh

echo $(date)

input_dir=/input
output_dir=/output

cwas permutation_test -i $output_dir/categorization/CWAS_input.categorization_result.zarr \
-o_dir $output_dir/permutation_test/ \
-s $input_dir/phenotype_table.tsv \
-p 35 \
-n 1000 \
-b --use_n_carrier

echo $(date)
echo "Done"
