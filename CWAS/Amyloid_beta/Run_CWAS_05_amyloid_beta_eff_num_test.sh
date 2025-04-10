#!/bin/sh

echo $(date)

input_dir=/input
output_dir=/output

cwas effective_num_test \
-i $output_dir/correlation/CWAS_input.correlation_matrix.zarr \
-c_count $output_dir/CWAS_input.category_counts.txt \
-o_dir $output_dir/eff_num_test \
-s $input_dir/phenotype_table.tsv \
-ef

echo $(date)
echo "Done"

