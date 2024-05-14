#!/bin/sh
#SBATCH -J 00.preprocessing
#SBATCH -p cpu
#SBATCH -o /data1/mycho/WGS_AD_2011/log/QC-hail/batch1-7/%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=40:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END

module purge

echo "___________START___________" 
date

##01.concat & sort
echo "___________01.concat & sort___________" 
date
DIR=/data2/Oneomics/VCF_hg38.1-6/batch1-7/
prefix=SMC_WGS.batch1-7

fn_input=$DIR/chr/chr*.annotated.vcf.gz
fn_output=$DIR/${prefix}.vcf.gz

#bcftools concat $fn_input -Oz -a --threads 16 -o $fn_output
echo "__reheader__"
##https://github.com/samtools/bcftools/issues/1041
#bcftools view -h $DIR/${prefix}.vcf.gz > $DIR/${prefix}.header


#fn_input=$DIR/${prefix}.vcf.gz
#fn_output=$DIR/${prefix}.sorted.vcf.gz

#bcftools sort -Oz $fn_input -o $fn_output

##02.rename (WGS_0285 to WGS_0930, WGS_0693 to WGS_0929)
echo "___________02.rename___________" 
date

fn_sampleID=$DIR/reheader.list
fn_input=$DIR/${prefix}.vcf.gz
fn_output=$DIR/${prefix}.rename.vcf

#fn_input=$DIR/${prefix}.sorted.vcf.gz
#fn_output=$DIR/${prefix}.sorted.rename.vcf
echo "WGS_0285 WGS_0930" > $fn_sampleID
echo "WGS_0693 WGS_0929" >> $fn_sampleID

bcftools reheader -s $fn_sampleID  -o $fn_output $fn_input 

##03. remove sample (WGS0345)
echo "___________03. remove sample_WGS0345___________" 
date
fn_input=$DIR/${prefix}.rename.vcf
fn_output=$DIR/${prefix}.rename.rm_sample.vcf.gz

#fn_input=$DIR/${prefix}.sorted.rename.vcf
#fn_output=$DIR/${prefix}.sorted.rename.rm_sample.vcf.gz
vcftools --remove-indv "WGS0345" --gzvcf $fn_input --recode --recode-INFO-all --out $fn_output

#fn_input=$DIR/${prefix}.sorted.rename.rm_sample.vcf.gz
#fn_output=$DIR/${prefix}.sorted.rename.rm_sample.vcf.bgz
fn_input=$DIR/${prefix}.rename.rm_sample.vcf.gz
fn_output=$DIR/${prefix}.rename.rm_sample.vcf.bgz
gunzip -c $fn_input |bgzip > $fn_output
