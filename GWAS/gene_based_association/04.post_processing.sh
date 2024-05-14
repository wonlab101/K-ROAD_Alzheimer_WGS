#!/bin/bash
DIR=/data1/mycho/WGS_AD_2011/
fn="SMC_batch1-7_hg38.QCed.final_230808"
prefix=231220_binary_geneburden_phenotype_update_240306


inputDIR=/data2/tmp_backup/minyoung/$prefix/step2/
#inputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/step2/
outputDIR=/data2/tmp_backup/minyoung/$prefix/merge/
#outputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/merge/
bfile=/data2/Oneomics/WGS/batch1-7/SMC_batch1-7_hg38.QCed.final_230808

phenolist=" DX_DAT_CU        RCL40_visual"

# make output PLOT directory
if [ ! -d $outputDIR/plot/ ]
then
         mkdir -p $outputDIR/plot/
 fi
 # 01. merge output

echo "____START merge___"
for pheno in $phenolist
do

for chr in {1..22}
do

        cat $inputDIR/step2.chr${chr}.gene_burden_${pheno}.regenie \
                | grep -v "^#" \
                | sed 's/ /\t/g' \
                >> $outputDIR/step2.merged.gene_burden.${pheno}.regenie

done

echo "____START log10Pval to PVAL___"
input=$outputDIR/step2.merged.gene_burden.${pheno}.regenie
Rscript pipe_log10Pval2Pval.R $input ${input}.pval

done


#02. Manhattan plot


echo "____START Manhattan plot___"

 for pheno in $phenolist
 do
         input=$outputDIR/step2.merged.gene_burden.${pheno}.regenie
         output=$outputDIR/plot/""$pheno".Manhattanplot"

         Rscript /data/kbs/scripts/ManhattanGG/ManhattanGG.R --assoc ${input}.pval --y-lim 7 \
                 --chr CHROM --pos GENPOS --snp ID --pval PVAL \
                 --genomewideline 2.9e-06 \
                 --suggestiveline 1e-05 --highlight-color yes \
                 --out ${output}.ylim7


 done
#03. QQ only

echo "____START QQ plot___"
 for pheno in $phenolist
 do
         input=$outputDIR/step2.merged.gene_burden.${pheno}.regenie
         output=$outputDIR/plot/""$pheno".qq"
        
         echo "##############################"$phenoname"____QQ"
         Rscript /data/kbs/scripts/ManhattanGG/ManhattanGG.R \
                 --assoc ${input}.pval --qq-only ci --snp ID --pval PVAL  --out $output

 done


 #04. significant

echo "____START significant___"

 for pheno in $phenolist
 do
         input_sig=$outputDIR/step2.merged.gene_burden.${pheno}.regenie.pval
         output_sig=$outputDIR/step2.merged.gene_burden.${pheno}.regenie.pval.1e02
        
         echo "##############################"$phenoname"____sig"
         
         head -n 1 $input_sig > $output_sig
         awk '{if($14<1e-02) print}' $input_sig >> $output_sig
        
 done





