#!/bin/sh

echo $(date)

echo "func_annot_func_score_gset.with_all"

input_dir=/data2/cwas_paper/cwas_input/SMC_input
sample_list=$input_dir/pheno_SMC_batch1-7_DXnoMCI_969_240422.tsv
out_dir=/data2/cwas_paper/cwas_output/SMC_output_969.annotation_v6.1

cat_info_path=$out_dir/feature_selection/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.train_prop80.rare_cutoff2.seed42.func_annot_func_score_gset.with_all.positive_r2.20240423.category_info.txt
total_cat_res_path=$out_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.sorted.categorization_result.zarr
filtered_cat_save_name_path=$out_dir/categorization/SMC_batch1-7_RCL40_visual_rare_het_AF_DB_MAF0.1_969_20240422_annot6.1.train_prop80.rare_cutoff2.seed42.func_annot_func_score_gset.with_all.positive_r2.20240423.categorization_result.zarr

cwas risk_score -i $filtered_cat_save_name_path \
-o_dir $out_dir/risk_score_after_fs \
-s $sample_list \
-c $cat_info_path \
--domain_list noncoding,intergenic,promoter,splice,intron,3primeutr,5primeutr,coding,ptv,missense \
--use_n_carrier \
-thr 2 \
-p 50 \
-tf 0.8 \
-t tf0.8 \
-pt 5,5

echo $(date)
echo "Done"