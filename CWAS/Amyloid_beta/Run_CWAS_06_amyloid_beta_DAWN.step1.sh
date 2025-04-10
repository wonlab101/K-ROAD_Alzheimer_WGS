#!/bin/sh

echo $(date)

output_dir=/output
cat_info_path=$output_dir/CWAS_input.category_info.txt 
cor_mat=$output_dir/correlation/CWAS_input.correlation_matrix.zarr 

counts=$output_dir/CWAS_input.category_counts.txt 
perm_test=$out_dir/permutation_test/CWAS_input.permutation_test.txt.gz

cwas effective_num_test -i $cor_mat \
-o_dir /output/eff_num_test \
-thr 10 \
-if corr -n 10000 \
--domain_list intergenic \
-t ef \
--category_set $cat_info_path \
-c_count $counts

echo $(date)
echo "Done"
