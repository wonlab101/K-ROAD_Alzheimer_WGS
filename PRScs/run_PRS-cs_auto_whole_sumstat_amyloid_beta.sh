#!/bin/sh
#SBATCH -J PRScs
#SBATCH -p cpu
#SBATCH -o log/log_%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=40:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END

ml purge
ml wonlab
ml ohpc

N_THREADS=10

export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

prefix=Carlos_NHW
dir_PRSCS=/data1/software/PRScs/
f_bed_orig=/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/SMC_batch1-7_hg38.QCed.final_230808.rsid_annovar
dir_output=/data1/mycho//WGS_AD_2011/10.PRScs/Calros/UKBB_all_snp/
f_bed=/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/data/SMC_batch1-7_hg38.QCed.final_230808.rsid_annovar.rsonly
f_rsid=/data1/mycho/GWAS_AD_protective/10.PRS/NG_PRS_ID_rsID.txt

f_NG_gwas=/data1/mycho/sumstat/Alzheimer/carlos/AmyloidPET_NHW_cohorts14_n11816_hg38_meta_WashU.txt
f_rsid_sst=/data1/mycho/sumstat/Alzheimer/carlos/AmyloidPET_NHW_cohorts14_n11816_hg38_meta_WashU.txt.anno
f_sst=/data1/mycho/WGS_AD_2011/10.PRScs/Calros/UKBB_all_snp/calros_AmyloidPET_NHW_cohorts14_n11816_hg38_meta_WashU.all.txt

##00.make sst file
#add rsID_ANNOVAR
cat header > $f_sst
cat $f_NG_gwas".anno"  | awk -F'\t' '$1 != 19 || ($4 < 43895848 || $4 > 45996742)' | grep "rs" | grep -v "CHR" | awk 'BEGIN{OFS="\t"} {print $15, $4, $5, $6,$8}' >> $f_sst


##01. PRScs
python $dir_PRSCS/PRScs.py \
        --ref_dir=$dir_PRSCS/ldblk/ldblk_ukbb/ldblk_ukbb_eur \
        --bim_prefix=$f_bed \
        --sst_file=$f_sst \
        --n_gwas=11816 \
        --phi=1e-2 \
        --out_dir=$dir_output
#IGAP n_gwas 63926
#NG n_GWAS 487511

cat $dir_output/*pst*txt > $dir_output/score.txt

##02. SCORE

# extract
cat $dir_output/score.txt | cut -f2 > $dir_output/$prefix".snplist"

plink --bfile $f_bed \
        --extract $dir_output/$prefix".snplist" \
        --make-bed --out $dir_output/$prefix

##PRScs output - hg19, plink - hg38
python /data/kbs/scripts/gwas/mk_bim2abc.py $dir_output/$prefix".bim" 
plink --threads 4 --bed $dir_output/$prefix".bed" --bim $dir_output/$prefix".abc.bim" --fam $dir_output/$prefix".fam" --make-bed --out $dir_output/$prefix".id"

python pipe_score2abc_hg38.py $dir_output/score.txt $dir_output/$prefix".bim" $dir_output/score.abc.txt

plink --bfile $dir_output/$prefix".id" \
        --score $dir_output/score.abc.txt 2 4 6 sum \
        --out $dir_output/PRS_cs.final.score_sum
