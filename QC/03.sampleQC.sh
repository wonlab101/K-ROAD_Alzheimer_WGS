#!/bin/sh
#SBATCH -J 03.sampleQC
#SBATCH -p cpu
#SBATCH -o /data1/mycho/WGS_AD_2011/log/QC-hail/batch1-7/%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=60
#SBATCH --time=60:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END

# Retrieve the job ID associated with the job name 02.genotypeQC
#JOBID=$(squeue --name=02.genotypeQC --format=%A --noheader)

# Set the job dependency using the obtained JOBID
####XXX SBATCH --dependency=afterok:${JOBID}

module purge
ml wonlab
ml ohpc

CONDA_PATH=/data1/software/anaconda3
ENV_NAME=hail

ENV_PATH=$CONDA_PATH/envs/$ENV_NAME 
source $CONDA_PATH/bin/activate $ENV_PATH 

echo "___________START___________"
date

###01. Run hail
python pipe_03.sampleQC.py

###02. KING
DIR=/data1/mycho/WGS_AD_2011/4.VCF.QCed.set1-7/
fn=SMC_batch1-7_hg38.sampleQC2_AT5.VQ1_230705
ln=SMC_batch1-7_hg38.sampleQC2_AT5.VQ1_230705
:<<"END"
# make log directory
if [ ! -d $DIR/log/ ]
then
                mkdir -p $DIR/log/
fi

if [ ! -d $DIR/output/KING/ ]
then
                mkdir -p $DIR/output/KING/
fi

if [ ! -d $DIR/output/PCA/ ]
then
                mkdir -p $DIR/output/PCA/
fi

### 0. make bed/bim/fam
echo "********make bed/bim/fam************"
date

plink2 --vcf $DIR/data/${fn}.vcf.bgz --threads 60 --const-fid  --make-bed --out $DIR/data/${fn}


### 1. KING relatedness
echo "**************KING start************"
date

plink --threads 60 --bfile $DIR/data/${fn} --chr 1-22 --maf 0.05 --geno 0.01 --make-bed --out $DIR/output/KING/${fn}.king > $DIR/log/${ln}_KING_1_plink.log

##KING options
## --cpus : specifies the number of CPU cores to be used for parallel computing
## --unrelated : a handy option to extract a list of unrelated individuals.
## --degree : specifies the degree of relatedness _an int number_, to go with KING analysis.

king -b $DIR/output/KING/${fn}.king.bed --rplot --cpus 60 --unrelated --degree 2 --prefix $DIR/output/KING/${fn}.king > $DIR/log/${ln}_KING_2_king.log


king -b $DIR/output/KING/${fn}.king.bed --rplot --cpus 60 --related --degree 2 --prefix $DIR/output/KING/${fn}.king.related > $DIR/log/${ln}_KING_2_king.related.log
king -b $DIR/output/KING/${fn}.king.bed --rplot --cpus 60 --duplicate --prefix $DIR/output/KING/${fn}.king.duplicate > $DIR/log/${ln}_KING_2_king.duplcate.log

plink --threads 60 --bfile $DIR/data/${fn} --remove $DIR/output/KING/${fn}.kingunrelated_toberemoved.txt --make-bed --out $DIR/output/${fn}.unrelated > $DIR/log/${ln}_KING_3_rm-unrelated.log

cat $DIR/output/KING/${fn}.king.related.kin | awk '$9 > 0.0884 {print}' > $DIR/output/KING/${fn}.king.related.kin.related

echo "***********KING END************"
date

### 2. PCA relatedness
echo "**************PCA start************"
date


# 2.1. Pruning

awk '{OFS="\t"}{print $1,$2,$3,$4,$5,"1"}' $DIR/output/${fn}.unrelated.fam > $DIR/output/${fn}.unrelated.aff1.fam


plink --bed $DIR/output/${fn}.unrelated.bed --bim $DIR/output/${fn}.unrelated.bim --fam $DIR/output/${fn}.unrelated.aff1.fam --chr 1-22 --maf 0.05 --geno 0.02 --make-bed --out $DIR/output/PCA/${fn}.unrelated.pca \
                >  $DIR/log/${ln}_PCA_1_replace_famfile.log

plink --bfile $DIR/output/PCA/${fn}.unrelated.pca --indep-pairwise 50 5 0.5 --out $DIR/output/PCA/${fn}.unrelated.pca \
                >  $DIR/log/${ln}_PCA_2_pruning.log

plink --bfile $DIR/output/PCA/${fn}.unrelated.pca --extract $DIR/output/PCA/${fn}.unrelated.pca.prune.in --make-bed --out $DIR/output/PCA/${fn}.unrelated.pca.pruned \
                >  $DIR/log/${ln}_PCA_3_makebed.log

echo "**************PCA_____END_____1.Pruning************"
# 2.2. pca

plink2 --bfile $DIR/output/PCA/${fn}.unrelated.pca.pruned --freq --pca approx var-wts --threads 60 --out $DIR/output/PCA/${fn}.unrelated.pca \
                >  $DIR/log/${ln}_PCA_4.runPCA.log


echo "**************PCA_____END_____2.PCA************"

END
#2.3. pca plot


Rscript /data1/mycho/WGS_AD_2011/script/QC-hail/KING_PCA/pipe_pca_plotting.R $DIR/output/PCA/${fn}.unrelated.pca

echo "***********PCA END************"
date
