#!/bin/bash
DIR=/data1/mycho/WGS_AD_2011/
fn="SMC_batch1-7_hg38.QCed.final_230808"
prefix=231220_binary

inputDIR1=$DIR/5.Association.1824.hg38/regenie/plink_chr_230817/
inputDIR2=$DIR/5.Association.1824.hg38/regenie/$prefix/step1/
outputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/step2/
cov=$DIR/1.data/phenotype/WGS_1824_phenotype_update.231218.txt

if [ ! -d $outputDIR ]
then
       mkdir -p $outputDIR
fi


for chr in {1..22}
do

        echo "*************START CHR" $chr


regenie_v2.2.4_hpc \
        --step 2 \
        --bed $inputDIR1/${fn}.chr$chr \
        --extract $inputDIR1/snp_pass_MAC31.snplist \
        --covarFile $cov \
        --covarCol age,PC{1:10} \
        --catCovarList sex,batch \
        --phenoFile $cov \
        --phenoCol DX_DAT_CU,RCL40_visual \
        --threads 40 \
        --bsize 100 \
        --bt \
        --firth --approx \
        --pThresh 0.01 \
        --pred $inputDIR2/step1.WGS_1824.hg38.chr${chr}_pred.list \
        --out $outputDIR/step2.WGS_1824.hg38.chr${chr}.single_variants


done


