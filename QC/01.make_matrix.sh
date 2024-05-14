#!/bin/sh
#SBATCH -J 01.make_matrix
#SBATCH -p cpu
#SBATCH -o /data1/mycho/WGS_AD_2011/log/QC-hail/batch1-7/%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=40:00:00
#SBATCH --mail-user=minyoungcho93@gmail.com
#SBATCH --mail-type=BEGIN,FAIL,END


module purge
CONDA_PATH=/data1/software/anaconda3
ENV_NAME=hail

ENV_PATH=$CONDA_PATH/envs/$ENV_NAME 
source $CONDA_PATH/bin/activate $ENV_PATH 

echo "___________START___________"
date


python pipe_01.make.matrix.py 
