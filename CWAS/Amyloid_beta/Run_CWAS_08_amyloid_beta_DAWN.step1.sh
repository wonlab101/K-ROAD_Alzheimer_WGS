#!/bin/sh

echo $(date)

out_dir=/data1/cwas_output/SMC_output.annotation_v6.1
cat_info_path=$out_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_info.txt
cor_mat=$out_dir/correlation/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.correlation_matrix.all.zarr

counts=$out_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.category_counts.txt
perm_test=$out_dir/permutation_test/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_1559_231220_annot6.0.train_prop80.rare_cutoff3.seed42.func_annot_func_score_gset.with_all.positive_r2.20240416.permutation_test.txt.gz

cwas effective_num_test -i $cor_mat \
-o_dir $out_dir/eff_num_test \
-thr 10 \
-if corr -n 10000 \
--domain_list intergenic \
-t efcutoff10 \
--category_set $cat_info_path \
-c_count $counts

echo $(date)
echo "Done"
