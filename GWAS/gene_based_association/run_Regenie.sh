#!/bin/sh
#SBATCH -J regenie_geneburden_231220_240306
#SBATCH -p cpu
#SBATCH -o /data1/mycho/WGS_AD_2011/log/association/1824.hg38/regenie//%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END


ml purge
ml wonlab
ml ohpc


echo "start 03.step2.regression"
date
sh 03.step2.regression.sh 


echo "start 04.post_processing"
date
sh 04.post_processing.sh 

echo "___FINISHED___"

