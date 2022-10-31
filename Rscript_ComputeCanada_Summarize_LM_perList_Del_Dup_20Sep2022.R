
#-------------------------------------------------------------------------------
# Summary: R-script to summarize the regression models saved
#-------------------------------------------------------------------------------

#-----
# # Summary: 
#          -> R-script to get summary DF for the regression model ran 
#          
# Key steps: 
#       Loop over gene-lists and create summary DF
#       
#-----
# Code:
# Author: Kuldeep Kumar
# Date: 20 September, 2022
# Version: 5
#
#-----
# UPDATE log: 
# .... Keep log of updates .....
#
#-----
# Input parameters:
#       1) flag_cohort_C  : a variable that loads Cohort specific files for CNV calls and Phenotype -> could be other file descriptor
#       2) flag_phenotype_C  : name of the phenotype variable to be used in Regression
#       3) filename_GeneList_C : input gene-list name (RData file of gene-lists)
#       
#  ** ADD other parameters if needed
#
#-----
# Output (saved as RData files in results_dir):
#        1) DF with summary of Regression models across List items (one for DEL and one for DUP)  
#
#-----


#-----------------------------------------------------------------------------
# Directory and paths: SPECIFY your own paths + directory structure
#-----------------------------------------------------------------------------

#----- Project Directory and other paths ---------------------------------
flag_local = 2  # 1: Local system; Else: Compute Canada Scratch Space or project directory path
if( flag_local == 1){
  project_dir <- "D:/SSDN_KD/Temp_KD/Temp_CNV_GWAS_oct2022"
} else {
  project_dir <- "/scratch/kkumar/Analysis/Temp_CNV_GWAS_oct2022"
}

code_dir <- paste0(project_dir,"/code")
data_dir <- paste0(project_dir,"/data")
results_dir <- paste0(project_dir,"/results")

#-----  Set Working directory (Code directory)
in_wd <- code_dir
setwd(in_wd)

#-----------------------------------------------------------------------------
# Packages
#-----------------------------------------------------------------------------

library(dplyr)
library(stringr)
library(R.utils)
#-----------------------------------------------------------------------------
# PART A: Parameters and variables
#-----------------------------------------------------------------------------

# Sample SLURM Job array submission:
# sbatch --time=30:00 script_SLURM_call_Rscript_Summarize_perListLM_oct2022.sh UKBB ZScore_IQ_adj_age_sex RData_file_lists_unique_gene_per_ROI_WholeGenome_15kAHBA 

# example assigned variables
#flag_cohort_C <- "UKBB" #"Mega"
#flag_phenotype_C <- "ZScore_IQ_adj_age_sex"
#filename_GeneList_C <- "RData_file_lists_unique_gene_per_ROI_WholeGenome_15kAHBA"


#-----------------------------------------------------------------------
#------ parse input command-line arguments

#------
# NOTE: we append _I, _N, _L, and _C to determine type of argument to be 
#             integer, numeric, logical, and string/character
# NOTE 2: the variables will then be assigned => global variables
#-----

cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

print(cli)

for (e in args) {
  argname <- e[1]
  if (! is.na(e[2])) {
    argval <- e[2]
    ## regular expression to delete initial \" and trailing \"
    argval <- gsub("(^\\\"|\\\"$)", "", argval)
  }
  else {
    # If arg specified without value, assume it is bool type and TRUE
    argval <- TRUE
  }
  
  print(e[1])
  print(e[2])
  
  # Infer type from last character of argname, cast val
  type <- substring(argname, nchar(argname), nchar(argname))
  if (type == "I") {
    argval <- as.integer(argval)
  }
  if (type == "N") {
    argval <- as.numeric(argval)
  }
  if (type == "L") {
    argval <- as.logical(argval)
  }
  if (type == "C") {
    #argval <- toString(dQuote(argval))
    argval <- toString(argval)
  }
  assign(argname, argval,env = .GlobalEnv)
  cat("Assigned", argname, "=", argval, "\n")
}

#-----------------------------------------------------------------------

#----- set experiment name => to be used when saving results as RData files
exp_name <- flag_cohort_C

#------ Option: ENSEMBL_ID or Gene_NAME
# Flag for ENSEMBL_ID or Gene-name => to be used for matching;
#      **SET this flag based on the input gene-list format
#      ** NOTE: gene_by_CNV file should have both columns: 
#               "gene_id" column for ENSEMBL_ID
#               "name_id" column for gene-Names

flag_ENSEMBL_GeneNAME <- 2
if(flag_ENSEMBL_GeneNAME == 1){
  gene_ID_colname <- "gene_id"
} else {
  gene_ID_colname <- "name_id"
}

#-----------------------------------------------------------------------------
# PART B: gene lists
#-----------------------------------------------------------------------------

#------ 4. READ / Load gene lists ==> Model would be run one list at a time
load(paste0(data_dir,"/",filename_GeneList_C,".RData"))   # **CHANGE the list name for other lists; 

#GeneList_filename <- paste0(data_dir,"/RData_file_lists_unique_gene_per_ROI_WholeGenome_15kAHBA.RData")  # Replace this with Any other list
#------ assign the loaded list to a new variable
list_Input_GeneSets <- list_unique_gene_per_ROI_WholeGenome   

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

#-----------------------------------------------------------
# function: Summarize regression output per list to a data-frame

fSummarize_LM_perList <- function(results_dir,exp_name,input_TYPE,n_gene_lists){
  
  #-------------------
  # Function Summary:
  #     -> summary of per-List regression models
  #
  #-------------------
  
  list_index <- c(1:n_gene_lists)
  Estimate <- c()
  StdError <- c()
  tValue <- c()
  pValue <- c()
  AdjRsquared <- c()
  Fstat <- c()
  DF <- c()
  
  for( select_index_I in c(1:n_gene_lists)){
    
    model_RData_filename <- paste0(results_dir,"/model_",input_TYPE,"_",exp_name,"_listIndex_",select_index_I,".RData")
    
    load(model_RData_filename)
    
    if(input_TYPE == "DEL"){
      lm_summary <- summary(model_DEL)
    } else {
      lm_summary <- summary(model_DUP)
    }
  
    #------Append summary stats
    Estimate <- c(Estimate,lm_summary$coefficients[2,1])
    StdError <- c(StdError,lm_summary$coefficients[2,2])
    tValue <- c(tValue,lm_summary$coefficients[2,3])
    pValue <- c(pValue,lm_summary$coefficients[2,4])
    
    AdjRsquared <- c(AdjRsquared,lm_summary$adj.r.squared)
    Fstat <- c(Fstat,lm_summary$fstatistic[[1]])
    DF <- c(DF,lm_summary$df[[2]])
  }
  
  # create a summary df
  df_summary_lm_perList <- data.frame(list_index = list_index,
                                      Estimate = Estimate,
                                      StdError = StdError,
                                      tValue = tValue,
                                      pValue = pValue,
                                      AdjRsquared = AdjRsquared,
                                      Fstat = Fstat,
                                      DF = DF)
  
  # return the summary df
  return(df_summary_lm_perList)
  
}



#-----------------------------------------------------------------------------
# Summarize Regression models perList for DEL and DUP
#-----------------------------------------------------------------------------

n_gene_lists <- length(list_Input_GeneSets)

#----- 1. Regression models for DEL
input_TYPE = "DEL"
df_summary_lm_perList_DEL <- fSummarize_LM_perList(results_dir,exp_name,input_TYPE,n_gene_lists)

colnames(df_summary_lm_perList_DEL)[c(2:ncol(df_summary_lm_perList_DEL))] <- paste0(colnames(df_summary_lm_perList_DEL)[c(2:ncol(df_summary_lm_perList_DEL))],"_DEL")

#----- 2. Regression models for DUP
input_TYPE = "DUP"
df_summary_lm_perList_DUP <- fSummarize_LM_perList(results_dir,exp_name,input_TYPE,n_gene_lists)

colnames(df_summary_lm_perList_DUP)[c(2:ncol(df_summary_lm_perList_DUP))] <- paste0(colnames(df_summary_lm_perList_DUP)[c(2:ncol(df_summary_lm_perList_DUP))],"_DUP")

#------ merge df: DEL and DUP
df_summary_lm_perList_DEL_DUP <- merge(df_summary_lm_perList_DEL,df_summary_lm_perList_DUP,by = "list_index")

# ----- save Model _DUP as RData file
save(df_summary_lm_perList_DEL_DUP,file = paste0(results_dir,"/df_summary_lm_perList_DEL_DUP_",exp_name,".RData"))
write.csv(df_summary_lm_perList_DEL_DUP,file = paste0(results_dir,"/df_summary_lm_perList_DEL_DUP_",exp_name,".csv"),row.names = FALSE)


#-----------------------------------------------------------------------------
# END
#-----------------------------------------------------------------------------
#

