#!/bin/bash
#SBATCH -J temp_meta
#SBATCH -p general
#SBATCH -A r00574
#SBATCH -o log_%x_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=90:00
#SBATCH --mail-user=minycho@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END

module load python
module load r

fn_1=/N/u/minycho/Quartz/slate_my/Projects/temp/meta/chip/step2.WGS_1824.hg38.merged.single_variants.DX_CU_DAT.regenie.pval

fn_2=/N/u/minycho/Quartz/slate_my/Projects/temp/meta/wgs/step2.WGS_1824.hg38.merged.single_variants.DX_DAT_CU.regenie.pval

fn_3=/N/u/minycho/Quartz/slate_my/Projects/temp/meta/wgs/NCGG_AD_WGS_GWAS.sumstat.MAC90.hg38.txt

output=/N/u/minycho/Quartz/slate_my/Projects/temp/meta/output//WGS1824_SMCGWAS_japanese_DX_CU_DAT.231228.hg38

python pipe_extract_commmon.py $fn_1".metalinput" $fn_2".metalinput" $fn_3".metalinput" $output".SMCCHIP.common.txt" $output".SMCWGS.common.txt" $output".JapaneseWGS.common.txt"

##metal
/N/u/minycho/Quartz/slate_my/tools/metal/generic-metal/metal pipe_metal.1.cmd
## Manhattan plot
input=$output".1.tbl"

python pipe_metal_after_formatting.py $input
Rscript ManhattanGG.R \
        --assoc $input".formatting" \
        --y-lim 10 \
        --chr CHR --pos POS --snp MarkerName --pval P-value \
        --suggestiveline 1e-06 --highlight-color yes \
        --out ${output}.ylim10.1e06



##significant
head $output".1.tbl.formatting" -n1 > $output".1.tbl.formatting.1e05"
awk '{if($6<1e-05) print}' $output".1.tbl.formatting" >> $output".1.tbl.formatting.1e05"

head $output".1.tbl.formatting" -n1 > $output".1.tbl.formatting.1e06"
awk '{if($6<1e-06) print}' $output".1.tbl.formatting" >> $output".1.tbl.formatting.1e06"
