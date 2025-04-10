#!/bin/sh

echo $(date)

input_dir=/input
output_dir=/output

cwas burden_shift \
-i $output_dir/CWAS_input.burden_test.txt \
-b $output_dir/permutation_test/CWAS_input.binom_pvals.txt.gz \
-c_info $output_dir/CWAS_input.category_info.txt \
-c_count $output_dir/CWAS_input.category_counts.txt \
-o_dir $output_dir/burden_shift/ \
-c_cutoff 9

echo $(date)
echo "Done"
