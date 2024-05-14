#!/bin/bash
prefix=231220_binary

DIR=/data1/mycho/WGS_AD_2011/
inputDIR=$DIR/5.Association.1824.hg38/regenie/plink_chr_230817/
outputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/step1
cov=$DIR/1.data/phenotype/WGS_1824_phenotype_update.231218.txt


fn="SMC_batch1-7_hg38.QCed.final_230808"

if [ ! -d $outputDIR ]
then
       mkdir -p $outputDIR
fi


for chr in {1..22}
do

        echo "*************START CHR" $chr


regenie_v2.2.4_hpc \
    --step 1 \
    --bed $inputDIR/${fn}.chr$chr \
    --extract $inputDIR/snp_pass_MAC31.snplist \
    --covarFile $cov \
    --covarCol age,PC{1:10} \
    --catCovarList sex,batch \
    --phenoFile $cov \
    --phenoCol DX_DAT_CU,RCL40_visual \
    --bsize 100 \
    --bt --lowmem \
    --niter 500 \
    --loocv \
    --lowmem-prefix $DIR/temp/tmp_rg_binary_cov_m1_230501 \
    --threads 40 \
    --out $outputDIR/step1.WGS_1824.hg38.chr$chr

done
