#!/bin/sh

echo $(date)

out_dir=/output
cat_info_path=$out_dir/CWAS_input.category_info.txt 
cor_mat=$out_dir/correlation/CWAS_input.correlation_matrix.zarr 
counts=$out_dir/CWAS_input.category_counts.txt 
perm_test=$out_dir/permutation_test/CWAS_input.permutation_test.txt.gz

lambdaval=3
corrval=0.12
countval=10
seed=57

cwas dawn \
-e $out_dir/eff_num_test/CWAS_input.intergenic.eig_vecs.ef.zarr \
-c $cor_mat \
-P $perm_test \
-o_dir $out_dir/dawn/ \
-k 194 \
-s $seed \
-T exact \
-t intergenic.C${countval}_R${corrval}_S2_L${lambdaval}_seed${seed} \
-c_count $counts \
-C $countval \
-R $corrval \
-S 2 \
-p 36 \
--lambda $lambdaval

echo $(date)
echo "Done"
