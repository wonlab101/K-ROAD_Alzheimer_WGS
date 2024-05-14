#!/bin/bash
DIR=/data1/mycho/WGS_AD_2011/
fn="SMC_batch1-7_hg38.QCed.final_230808"
prefix=231220_binary_geneburden_phenotype_update_240306


inputDIR1=$DIR/5.Association.1824.hg38/regenie/plink_chr_230817/
inputDIR2=$DIR/5.Association.1824.hg38/regenie/231220_binary/step1/
outputDIR=/data2/tmp_backup/minyoung/$prefix/step2
#outputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/step2/
cov=$DIR/1.data/phenotype/WGS_1824_phenotype_update.231218.txt


if [ ! -d $outputDIR ]
then
       mkdir -p $outputDIR
fi

annopath=/data2/Oneomics/WGS/batch1-7/VEP/gene_burden/SMC1824_hg38/
anno=$annopath/SMC_batch1-7_hg38.QCed.final_230808.info.merged.vep_table.anno.rmdup
setlist=$annopath/SMC_batch1-7_hg38.QCed.final_230808.info.merged.vep_table.setlist
mask=$annopath/SMC_batch1-7_hg38.QCed.final_230808.mask


for chr in {1..22}
do

        echo "*************START CHR" $chr


regenie_v2.2.4_hpc \
        --step 2 \
        --bed $inputDIR1/${fn}.chr$chr \
        --covarFile $cov \
        --covarCol age,PC{1:10} \
        --catCovarList sex,batch \
        --phenoFile $cov \
        --phenoCol DX_DAT_CU,RCL40_visual \
        --threads 10 \
        --bsize 100 \
        --bt \
        --firth --approx \
        --pThresh 0.01 \
        --pred $inputDIR2/step1.WGS_1824.hg38.chr${chr}_pred.list \
        --anno-file $anno \
        --set-list $setlist \
        --mask-def $mask \
        --write-mask \
        --aaf-bins 0.01 \
        --lowmem-prefix $DIR/temp/tmp_rg.maf01 \
        --out $outputDIR/step2.chr${chr}.gene_burden 

done
