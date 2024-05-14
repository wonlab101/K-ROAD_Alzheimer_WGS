#!/bin/bash
DIR=/data1/mycho/WGS_AD_2011/

prefix=231220_binary
inputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/step2/
outputDIR=$DIR/5.Association.1824.hg38/regenie/$prefix/merge/
bfile=$DIR/5.Association.1824.hg38/regenie/plink_chr_230817/SMC_batch1-7_hg38.QCed.final_230808
anno_input=/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/anno//SMC_batch1-7_hg38.QCed.final_230808

phenolist=" DX_DAT_CU     RCL40_visual"

# make output PLOT directory
if [ ! -d $outputDIR/plot/ ]
then
         mkdir -p $outputDIR/plot/
 fi

 # 01. merge output
echo "____START merge___"
for pheno in $phenolist
do

cat $inputDIR/step2.WGS_1824.hg38.chr1.single_variants_${pheno}.regenie \
        grep CHROM > $outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie

for chr in {1..22}
do

        cat $inputDIR/step2.WGS_1824.hg38.chr${chr}.single_variants_${pheno}.regenie \
                | grep -v CHROM \
                | sed 's/ /\t/g' \
                >> $outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie

done

echo "____START log10Pval to PVAL___"
input=$outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie
Rscript pipe_log10Pval2Pval.R $input ${input}.pval

done
#02. Manhattan plot


echo "____START Manhattan plot___"

 for pheno in $phenolist
 do
         input=$outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie
         output=$outputDIR/plot/WGS_1824.hg38.$pheno".Manhattanplot"
         Rscript pipe_log10Pval2Pval.R $input ${input}.pval
         
         Rscript /data/kbs/scripts/ManhattanGG/ManhattanGG.R --assoc ${input}.pval --y-lim 50 \
                 --chr CHROM --pos GENPOS --snp ID --pval PVAL \
                 --suggestiveline 1e-05 --highlight-color yes \
                 --out ${output}.ylim80


 done

#03. QQ only

echo "____START QQ plot___"
 for pheno in $phenolist
 do
         input=$outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie
         output=$outputDIR/plot/.WGS_1824.hg38.$pheno".qq"
        
         echo "##############################"$phenoname"____QQ"
         Rscript /data/kbs/scripts/ManhattanGG/ManhattanGG.R \
                 --assoc ${input}.pval --qq-only ci --snp ID --pval PVAL  --out $output

 done


 #04. significant

echo "____START significant___"

 for pheno in $phenolist
 do
         input_sig=$outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie.pval
         output_sig=$outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie.pval.1e05
        
         echo "##############################"$phenoname"____sig"
         
         head -n 1 $input_sig > $output_sig
         awk '{if($14<1e-05) print}' $input_sig >> $output_sig
        
 done

#05.lead SNP

echo "____START lead SNP___"
for pheno in $phenolist
do
        input_clump=$outputDIR/step2.WGS_1824.hg38.merged.single_variants.${pheno}.regenie.pval.1e05
        output_clump=${input_clump}.clump0.2

        #1) clumping
        plink --bfile $bfile --clump $input_clump --clump-snp-field ID --clump-field PVAL --clump-r2 0.2 --clump-kb 3000 --clump-p1 1e-05 --clump-p2 1e-05 --out $output_clump

        #2) make snp list
        snplist=${output_clump}.snplist
        awk 'NF {print $3}' ${output_clump}.clumped > $snplist

        #3) mark Lead
        output_lead=${input_clump}.leadSNP
        python3.6 /data1/mycho/script/GWAS/association_markLead.py ${input_clump} $snplist $output_lead


echo "____START ANNOVAR___"

        input=${input_clump}.leadSNP
        #sh annotation/01.make_frq_file.sh $input
        #sh annotation/02.annotation.sh $input
        

        annovar1=${anno_input}.anno.1.hg38_multianno.txt
        annovar2=${anno_input}.anno.2.hg38_multianno.txt

        output=${input}.anno

         python3.6 pipe_association_add_ANNOVAR_results.py $annovar1 $annovar2 $input $output


done


