#!/bin/sh
#SBATCH -J 01.preprocessing
#SBATCH -p cpu
#SBATCH -o /data1/mycho/WGS_AD_2011/log/association/1824.hg38/regenie/%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=40:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END

ml purge
ml wonlab
ml ohpc


DIR=/data1/mycho/WGS_AD_2011/

plink --bfile /data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/SMC_batch1-7_hg38.QCed.final_230808 \
        --threads 16 \
        --mac 31 \
        --keep /data1/mycho/WGS_AD_2011/1.data/phenotype/WGS_1824_phenotype_update_onlyADMCICU.230817_keeplist.txt \
        --write-snplist \
        --out $DIR/5.Association.1824.hg38/regenie/snp_pass_MAC31

for chr in {1..22}
do 
        plink --bfile /data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/SMC_batch1-7_hg38.QCed.final_230808 \
                --chr $chr --make-bed --threads 16 \
        --keep /data1/mycho/WGS_AD_2011/1.data/phenotype/WGS_1824_phenotype_update_onlyADMCICU.230817_keeplist.txt \
                --out /data1/mycho/WGS_AD_2011/5.Association.1824.hg38/regenie/plink_chr_230817/SMC_batch1-7_hg38.QCed.final_230808.chr$chr  

done
