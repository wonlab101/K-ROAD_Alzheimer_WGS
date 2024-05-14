#!/bin/sh
#SBATCH -J regenie_231220
#SBATCH -p cpu
#SBATCH -o /data1/mycho/WGS_AD_2011/log/association/1824.hg38/regenie//%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=100:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END

ml purge
ml wonlab
ml ohpc

CONDA_PATH=/data1/software/anaconda3        
ENV_NAME=REGENIE                              
ENV_PATH=$CONDA_PATH/envs/$ENV_NAME         
source $CONDA_PATH/bin/activate $ENV_PATH    

echo "start 01.data_processing"
date
echo "PASS___"
#sh 01.data_processing.sh > $logDIR/01.data_processing.log

echo "start 02.step1.fitting"
date
sh 02.step1.fitting_m1.sh 
echo "PASS__"

echo "start 03.step2.regression"
date
sh 03.step2.regression_m1.sh 


echo "start 04.post_processing"
date
sh 04.post_processing_m1.sh 

echo "___FINISHED___"



