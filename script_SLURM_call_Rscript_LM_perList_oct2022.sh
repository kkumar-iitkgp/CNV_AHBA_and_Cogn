#!/bin/bash

#SBATCH --account=rrg-jacquese
#SBATCH --job-name=CNV_GWAS
# #SBATCH --array=1-200
# #SBATCH --time=2:55:00
#SBATCH --mem=15000M
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#----------------------------------------------------------------------------------------------------
# directory
BASE_DIR=$SCRATCH/Analysis/Temp_CNV_GWAS_oct2022
cd ${BASE_DIR}/code


module load StdEnv/2020
module load gcc/9.3.0
module load r/4.1.0

export R_LIBS=~/.local/R/$EBVERSIONR/

#----------------------------------------------------------------------------------------------------
# read some parameter or variable values as input

INPUT_cohort_name=$1
INPUT_phenotype_name=$2
INPUT_GeneLIST_File=$3

echo "Starting task $SLURM_ARRAY_TASK_ID"

TASK_ID=${SLURM_ARRAY_TASK_ID}
#temp_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${INPUT_JOB_ARRAY_ARG_FILE})

#----------------------------------------------------------------------------------------------------
# 1. call Rscript one per List

R --vanilla -f Rscript_ComputeCanada_Run_CountSummary_SumScore_LM_perList_20Sep2022.R --args flag_cohort_C=${INPUT_cohort_name} flag_phenotype_C=${INPUT_phenotype_name} filename_GeneList_C=${INPUT_GeneLIST_File}  select_index_I=${TASK_ID} >log_Rscript_LM_perList_${SLURM_JOBID}.Rout 2>&1

# remove Rout log
#rm log_Step3_slurm_${SLURM_JOBID}.Rout

#----------------------------------------------------------------------------------------------------


