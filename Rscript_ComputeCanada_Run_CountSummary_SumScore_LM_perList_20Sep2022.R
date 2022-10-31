
#-------------------------------------------------------------------------------
# Summary: R-script to get sum-of-scores for an input gene-list and run regression model
#-------------------------------------------------------------------------------

#-----
# # Summary: 
#          -> R-script to get sum-of-scores and run regression model 
#          -> The script runs for one gene-list which is passed as "select_index_I" input from command line 
# Key steps: 
#       For a given gene-list and select_index_I run:
#       1) PART 1:OPTIONAL: Get count of genes (unique/total) observed in our data (Runs only once; for select_index_I = 1
#       2) PART 2: Using CNV-calls; phenotype info; and input gene-lists -> Get data frame with sum-of-scores + phenotype file 
#       3) PART 3: Run regression model -> save the model
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
#       4) select_index_I  : List index for which regression would be run (1 call per list); => PARALLEL Job submission
#  
#  ** ADD other parameters if needed
#
#-----
# Output (saved as RData files in results_dir):
#       1) Count summary: data-frame  with count of unique & total genes per gene-list
#       2) Sum-of_scores: data-frame with sum-of-score and phenotype info
#       3) Regression models: regression model (one for DEL and one for DUP)  
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
# sbatch --time=30:00 --array=1-180 script_SLURM_call_Rscript_LM_perList_oct2022.sh UKBB ZScore_IQ_adj_age_sex RData_file_lists_unique_gene_per_ROI_WholeGenome_15kAHBA 

# example assigned variables
#flag_cohort_C <- "UKBB" #"Mega"
#flag_phenotype_C <- "ZScore_IQ_adj_age_sex"
#filename_GeneList_C <- "RData_file_lists_unique_gene_per_ROI_WholeGenome_15kAHBA"
#select_index_I <- 1   => index of the list for which models will be run

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
# PART B: Read Data and gene lists
#-----------------------------------------------------------------------------

#------ **Suggestion for 1,2, & 3: save and then load as RData file => FASTER **-------#
# ----- *NOTE: these files were obtained from the base files from Guillaume after 
#              removing subjects based on filters suggested in Guillaume's email
#              + removing some columns that are not needed

#------ 1. CNV call info: gene by CNV file ----------------------------
###  Filtered_gene_by_CNV = read.csv(file=paste0(data_dir,"/Filtered_gene_by_CNV_Mega_UKBBMEGASPARKMSSNG_DEL_1_DUP_1_1nov2021.csv"))
Filtered_gene_by_CNV <- read.csv(paste0(data_dir,"/Filtered_gene_by_CNV_",flag_cohort_C,".csv"))

#------ 2. CNV call info: CNV by Individual file -----------------------------
### Filtered_CNV_by_IID <- read.csv(paste0(data_dir,"/CNVfiltered_overlapartifactFreq_UKBBMEGASPARKMSSNG_20211001.csv"))
Filtered_CNV_by_IID <- read.csv(paste0(data_dir,"/Filtered_CNV_by_IID_",flag_cohort_C,".csv"))

#------ 3. Phenotype file -----------------------
### Final_Phenotype_IND <- read.delim(paste0(data_dir,"/MEGA_MSSNG_SPARK_UKBB_pheno_scores_gene_100_20211001.tsv"))
# ------- remove subjects based on filters
Final_Phenotype_IND <- read.csv(paste0(data_dir,"/Final_Phenotype_IND_",flag_cohort_C,".csv"))


#------ 4. READ / Load gene lists ==> Model would be run one list at a time
load(paste0(data_dir,"/",filename_GeneList_C,".RData"))   # **CHANGE the list name for other lists; 

#GeneList_filename <- paste0(data_dir,"/RData_file_lists_unique_gene_per_ROI_WholeGenome_15kAHBA.RData")  # Replace this with Any other list
#------ assign the loaded list to a new variable
list_Input_GeneSets <- list_unique_gene_per_ROI_WholeGenome   

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

#-----------------------------------------------------------
# function 1: get Binary Data frame (for sum or count of score) for a gene-list and select index
fGeneList_to_BinaryDF_inputSelectIndex <- function(input_gene_list,select_index_I){
  
  #-------------------
  # Function Summary:
  # Aim: function get binary matrix with 3 columns -> 
  #       col 1: Gene names (or ENSEMBL IDs, depending on the list values)
  #       col 2: Binary column for all the genes -across-lists (all values will be 1)
  #       col 3: Binary column for the genes in select_index_I list => set genes within list to 1 else 0
  #-------------------
  
  # check: select_index_I <= length(input_gene_list)
  if( (select_index_I > length(input_gene_list)) | (select_index_I < 1 )){
    #----- Stop and throw error message: select_index_I can be a number between 1 & the length of the input list
    stop("Error: select index should be between 1 and the number of list elements!")
  } else {
    
    #----- Create a binary data-frame with information in all-genes-across-lists + genes-within-select-index-list
    #----- i. Vector of all the unique genes
    vec_unique_genes <- unique(unlist(input_gene_list))
    
    #----- ii. create df with gene-names; binary col for all genes; binary column for genes in select list (default value 0)
    df_gene_set_binary <- data.frame(gene_name=vec_unique_genes,
                                     Total_count = rep(1,length(vec_unique_genes)),
                                     Count_genes_within_list = rep(0,length(vec_unique_genes)))
    
    #----- iii. For genes in select index gene list set values of column  "bin_col_select_list" to 1
    df_gene_set_binary[df_gene_set_binary[,"gene_name"] %in% input_gene_list[[select_index_I]],"Count_genes_within_list"] <- 1
    
    #----- return input-gene-list converted to binary encoded data-frame
    return(df_gene_set_binary)
  }
  
}



#-----------------------------------------------------------
# # function 2: get Sum of scores or Count of genes within a gene-list for each individual =>
#               NOTE it's two level process: i. Sum of scores for genes within CNV and then ii. Sum of scores for CNVs within an individual
fSum_scores_or_genecount_per_Individual <- function(Filtered_CNV_by_IID,Filtered_gene_by_CNV,df_gene_score_or_category_list,input_TYPE,array_sum_score_names){
  
  #-------------------
  # Function Summary:
  # Aim: Get sum of scores per Individual using 
  #      i. CNV calls with ii. gene lists info (one list at a time)
  # NOTE: gene-list info is coded in "df_gene_score_or_category_list", which is a binary matrix & for a
  #        For a binary matrix the sum of scores is equivalent to count of genes.
  #-------------------
  
  #----- 1. Subset the CNV_by_Individual data-frame for input_TYPE (i.e. DEL or DUP)
  CNV_by_IID_subset = Filtered_CNV_by_IID[Filtered_CNV_by_IID["TYPE"] == input_TYPE,]
  
  #----- 2. Subset gene_by_CNV data-frame for input_TYPE (i.e. DEL or DUP)
  gene_by_CNV_subset = Filtered_gene_by_CNV[Filtered_gene_by_CNV[,"TYPE"] == input_TYPE,]
  
  #----- 3. gene_score_info: merge gene_by_CNV subset data-frame with the gene-list info ( => df_gene_score_or_category_list)
  # # # # ** NOTE: change "gene_name" with "ENSEMBL-ID" or other format for gene-IDs 
  # # # # *** CHECK: this has to be same format in both data frames
  gene_score_info = merge(gene_by_CNV_subset,df_gene_score_or_category_list,by.x=gene_ID_colname,by.y="gene_name")   # use "gene_id" for ENSEMBL_IDs
  
  #----- 4. Add new column "position" to "gene_score_info" based on CHR, START, STOP 
  #          => Will be used to make a sum of genes within a CNV
  gene_score_info[,"position"] <- paste0(gene_score_info[,"CHR"],"_",gene_score_info[,"START"],"_",gene_score_info[,"STOP"])
  
  # # # # Sum scores across all genes within a CNV -----------------------
  gene_score_info_select_cols <- gene_score_info[,c("position",array_sum_score_names)]
  SUM_GENES <- gene_score_info_select_cols %>% group_by(position) %>% summarise(across(everything(), list(sum)))
  

  #----- 5. Add new column "position" to "CNV_by_IID_subset" based on CHR, START, STOP 
  CNV_by_IID_subset[,"position"] <- paste0(CNV_by_IID_subset[,"CHR"],"_",CNV_by_IID_subset[,"START"],"_",CNV_by_IID_subset[,"STOP"])
  CNV_by_IID_subset_POS <- CNV_by_IID_subset[,c("individual","position")]
  
  
  #----- 6. Merge CNV_by_IID_pos and sum_genes
  POS_IND_score <- merge(CNV_by_IID_subset_POS,SUM_GENES,by="position",all=TRUE)
  
  # # # replace NA with 0 => count of 0 for sum across CNVs for an individual
  POS_IND_score[is.na(POS_IND_score)] <- 0
  
  #----- 7. Sum over all CNVs within an individual ---------------------------
  #  # # Sum across all CNVs within an Individual
  POS_IND_score_IID = POS_IND_score[, -which(names(POS_IND_score) %in% c("position"))]
  SUM_IND_SCORE = POS_IND_score_IID %>% group_by(individual) %>% summarise(across(everything(), list(sum)))
  
 
  # ----- Simplify column names ---------------------------
  colnames(SUM_IND_SCORE) <- c("individual",array_sum_score_names)
  
  #----- Return SUM of Scores data-frame => per-individual sum (first sum of genes within CNVs, and then Sum of CNVs within individual)
  return(SUM_IND_SCORE)
   
} 


#-----------------------------------------------------------
# # function 3: Data-frame for running regression model: Sum of Scores for each individual (DEL and DUP both merged) + Phenotype 
fPrepare_variables_for_regression_perList <- function(Filtered_CNV_by_IID,Filtered_gene_by_CNV,Final_Phenotype_IND,df_gene_score_or_category_list,array_sum_score_names){
  
  #-------------------
  # Function Summary:
  # Aim: function get / prepare variables for the regression models 
  #       => combine i. CNV calls with ii. gene lists info and iii. Phenotype 
  #       => get Sum of SCORES for DEL and DUP and then combine
  #       => Finally add "Phenotype" information (Which is an independent step)
  #
  #-------------NOTE: we replace "Sum-of-scores" NA values with 0 => they count as 0 score or 0 gene count
  #             However, we do not replace the NA values in phenotype column to 0; as this will add spurious phenotype points
  #-------------------
  
   
  #----- 1. DEL -> Sum across individual
  input_TYPE = "DEL"
  SUM_IND_SCORE_DEL <- fSum_scores_or_genecount_per_Individual(Filtered_CNV_by_IID,Filtered_gene_by_CNV,df_gene_score_or_category_list,input_TYPE,array_sum_score_names)
  # create column "Count_genes_not_in_list"
  SUM_IND_SCORE_DEL[,"Count_genes_not_in_list"] <- SUM_IND_SCORE_DEL[,array_sum_score_names[1]] - SUM_IND_SCORE_DEL[,array_sum_score_names[2]]
  names(SUM_IND_SCORE_DEL) <- c("individual",paste0(c(array_sum_score_names,"Count_genes_not_in_list"),"_",input_TYPE))
    
  #----- 2. DUP -> Sum across individual 
  input_TYPE = "DUP"
  SUM_IND_SCORE_DUP <- fSum_scores_or_genecount_per_Individual(Filtered_CNV_by_IID,Filtered_gene_by_CNV,df_gene_score_or_category_list,input_TYPE,array_sum_score_names)
  # create column "Count_genes_not_in_list"
  SUM_IND_SCORE_DUP[,"Count_genes_not_in_list"] <- SUM_IND_SCORE_DUP[,array_sum_score_names[1]] - SUM_IND_SCORE_DUP[,array_sum_score_names[2]]
  names(SUM_IND_SCORE_DUP) <- c("individual",paste0(c(array_sum_score_names,"Count_genes_not_in_list"),"_",input_TYPE))
  
  
  #----- 3. Merge SUM_IND_SCORE_DEL and SUM_IND_SCORE_DUP into a single DF
  POS_IND_score_merge <-  merge(SUM_IND_SCORE_DEL,SUM_IND_SCORE_DUP,by="individual",all=TRUE)
  
  ###----- replace NA values with 0
  POS_IND_score_merge[is.na(POS_IND_score_merge)] <- 0
  
  #----- 4. Add phenotype information
  FINAL_INDIVIDUAL_SCORES=merge(Final_Phenotype_IND,POS_IND_score_merge,by="individual",all=TRUE)
  
  # replace NA values with 0 for Count columns
  count_column_names <- c(paste0(c(array_sum_score_names,"Count_genes_not_in_list"),"_DEL"),paste0(c(array_sum_score_names,"Count_genes_not_in_list"),"_DUP"))
  for(temp_col in c(1:length(count_column_names))){
    FINAL_INDIVIDUAL_SCORES[is.na(FINAL_INDIVIDUAL_SCORES[,count_column_names[temp_col]]),count_column_names[temp_col]] <- 0
  }


  #----- Return the dataframe with phenotype + sum-of-scores for each individual 
  return(FINAL_INDIVIDUAL_SCORES)
  
} 


#-----------------------------------------------------------
# # function 4: Get Count of genes for each gene-list => genes observed in our data
#               NOTE it's two level process: i. take into account all genes and CNVs (including recurrence); ii) only individuals with phenotype
fCount_unique_totalRec_genes_across_genelist <- function(Filtered_CNV_by_IID,Filtered_gene_by_CNV,list_individuals_with_phenotype,list_Input_GeneSets,input_TYPE){
  
  #-------------------
  # Function Summary:
  # Aim: Count of unique / total-Reccurent genes from a gene-list observed in our data 
  #      i. CNV calls with ii. gene lists info (one list at a time)
  #
  # STEPS:
  #----- 1. Unique genes => genes within each list that are deleted or duplicated and are in gene_by_CNV file
  #----- 2. Total recurrent genes => count of all the genes that are deleted or duplicated => gene count after merging gene_by_CNV + CNV_by_Individual 
  #----- 3. Total individuals per gene-list => count of all individuals that are carry genes deleted or duplicated => individual count after merging gene_by_CNV + CNV_by_Individual 
  #----- 4.  N-genes_across_geneLists => all the genes considered in the experiment
  #
  #-------------------
  
  #----- 1. Subset the CNV_by_Individual data-frame for input_TYPE (i.e. DEL or DUP)
  CNV_by_IID_subset = Filtered_CNV_by_IID[Filtered_CNV_by_IID["TYPE"] == input_TYPE,]
  
  # # # Add new column "position" to "CNV_by_IID_subset" based on CHR, START, STOP 
  CNV_by_IID_subset[,"position"] <- paste0(CNV_by_IID_subset[,"CHR"],"_",CNV_by_IID_subset[,"START"],"_",CNV_by_IID_subset[,"STOP"])
  CNV_by_IID_subset_POS <- CNV_by_IID_subset[,c("individual","position")]
  
  
  #----- 2. Subset gene_by_CNV data-frame for input_TYPE (i.e. DEL or DUP)
  gene_by_CNV_subset = Filtered_gene_by_CNV[Filtered_gene_by_CNV[,"TYPE"] == input_TYPE,]
  
  #----- Add new column "position" to "gene_score_info" based on CHR, START, STOP 
  #          => Will be used to make a sum of genes within a CNV
  gene_by_CNV_subset[,"position"] <- paste0(gene_by_CNV_subset[,"CHR"],"_",gene_by_CNV_subset[,"START"],"_",gene_by_CNV_subset[,"STOP"])
  
  
  #----- 3. gene_by_CNV df for individuals with info
  All_gene_in_CNV_info = merge(CNV_by_IID_subset_POS,gene_by_CNV_subset,by="position",all=TRUE) 
  
  All_gene_in_CNV_info <- All_gene_in_CNV_info[complete.cases(All_gene_in_CNV_info[,c("individual","gene_id","name_id")]),]
  # Keep only those that have a phenotype
  All_gene_in_CNV_info <- All_gene_in_CNV_info[All_gene_in_CNV_info[,"individual"] %in% list_individuals_with_phenotype,]
  
  #  #  # NOTE: All_gene_in_CNV_info will be used for count
  
  
  #----- 4. Loop over gene-lists and get counts---------------------------
  n_gene_lists = length(list_Input_GeneSets)
  
  # Overall Experiment gene set: Count the number of unique genes in our experiment => Across lists
  vec_unique_genes <- unique(unlist(list_Input_GeneSets))
  array_Overall_count_GeneList_FullGeneSet <- rep(length(vec_unique_genes),n_gene_lists)
  
  # Overall data gene set: Count the number of genes observed in our data (overall set comparison => Across lists)
  vec_data_genes_observed <- All_gene_in_CNV_info[All_gene_in_CNV_info[,gene_ID_colname] %in% vec_unique_genes,gene_ID_colname]  # use "gene_id" for ENSEMBL_IDs
  array_Overall_count_GeneList_DataGenes_unique <- rep(length(unique(vec_data_genes_observed)),n_gene_lists)
  array_Overall_count_GeneList_DataGenes_Total <- rep(length(vec_data_genes_observed),n_gene_lists)
  array_Overall_count_GeneList_DataGenes_UniqueIndividual <- rep(length(All_gene_in_CNV_info[All_gene_in_CNV_info[,gene_ID_colname] %in% vec_unique_genes,"individual"]),n_gene_lists) # use "gene_id" for ENSEMBL_IDs
  
  # Per-List gene count
  array_count_perList_Input <- c()
  array_count_perList_unique <- c()
  array_count_perList_total <- c()
  array_count_perList_unique_individual <- c()
  for( loop_l in c(1:n_gene_lists)){
    
    # genes in a single list 
    vec_genes_temp_list <- unique(unlist(list_Input_GeneSets[[loop_l]]))
    array_count_perList_Input <- c(array_count_perList_Input,length(vec_genes_temp_list))
    
    # genes in the single list that are observed in our data
    vec_list_genes_observed <- All_gene_in_CNV_info[All_gene_in_CNV_info[,gene_ID_colname] %in% vec_genes_temp_list,gene_ID_colname]  # use "gene_id" for ENSEMBL_IDs
    array_count_perList_unique <- c(array_count_perList_unique, length(unique(vec_list_genes_observed)))
    array_count_perList_total <- c(array_count_perList_total, length(vec_list_genes_observed))
    array_count_perList_unique_individual <- c(array_count_perList_unique_individual, length(unique(All_gene_in_CNV_info[All_gene_in_CNV_info[,gene_ID_colname] %in% vec_genes_temp_list,"individual"]))) # use "gene_id" for ENSEMBL_IDs
    
  }
  
  # ----- Create data frame ---------------------------
  df_gene_count_per_list <- data.frame(Overall_count_GeneList_FullGeneSet = array_Overall_count_GeneList_FullGeneSet,
                                       Overall_count_GeneList_DataGenes_unique = array_Overall_count_GeneList_DataGenes_unique,
                                       Overall_count_GeneList_DataGenes_Total = array_Overall_count_GeneList_DataGenes_Total,
                                       Overall_count_GeneList_DataGenes_UniqueIndividual = array_Overall_count_GeneList_DataGenes_UniqueIndividual,
                                       count_perList_Input = array_count_perList_Input,
                                       count_perList_unique = array_count_perList_unique,
                                       count_perList_total = array_count_perList_total,
                                       count_perList_unique_individual = array_count_perList_unique_individual)
  
  #----- Return data-frame with count of genes
  return(df_gene_count_per_list)
  
} 



#-----------------------------------------------------------------------------
# PART 1: Count: Count genes observed in Data => unique / Total genes that are in data & within each list
#-----------------------------------------------------------------------------

#----- NOTE: This will be Run only ONCE => Here for select_index_I = 1, but summarizes acorss all the gene-list sets/items. 
#           * The single run covers info for all the gene-list sets. For example, for cortical regions it will have 180 values one per region
if( select_index_I == 1){
  
  list_individuals_with_phenotype <- Final_Phenotype_IND[!is.na(Final_Phenotype_IND[,flag_phenotype_C]),"individual"]
  
  #----- 1. DEL count summary
  input_TYPE <- "DEL"
  df_gene_count_per_list_DEL <- fCount_unique_totalRec_genes_across_genelist(Filtered_CNV_by_IID,Filtered_gene_by_CNV,list_individuals_with_phenotype,list_Input_GeneSets,input_TYPE)
  
  #----- 2. DUP count summary
  input_TYPE <- "DUP"
  df_gene_count_per_list_DUP <- fCount_unique_totalRec_genes_across_genelist(Filtered_CNV_by_IID,Filtered_gene_by_CNV,list_individuals_with_phenotype,list_Input_GeneSets,input_TYPE)
  
  #----- save RData file: Both DEL and DUP summary df in single RData file
  save(df_gene_count_per_list_DEL,df_gene_count_per_list_DUP,file = paste0(results_dir,"/summary_df_gene_count_per_list_DEL_DUP_",exp_name,".RData"))

}

#-----------------------------------------------------------------------------
# PART 2: Sum of Scores => prepare data-frame for running regression
#-----------------------------------------------------------------------------

#------ Binary data-frame with information on all the genes (Across all lists) + binary coding for the select_index_I list
df_gene_score_or_category_list <- fGeneList_to_BinaryDF_inputSelectIndex(list_Input_GeneSets,select_index_I)

#----- Name of the columns / scores to be summed 
#   NOTE: *MODIFY this for change of score names or ** For regression model 
array_sum_score_names <- c("Total_count","Count_genes_within_list")         # => c("bin_col_genes_all","bin_col_genes_select_list")

#----- Sum of Scores + Phenotype for each individual 
FINAL_INDIVIDUAL_SCORES <- fPrepare_variables_for_regression_perList(Filtered_CNV_by_IID,Filtered_gene_by_CNV,Final_Phenotype_IND,df_gene_score_or_category_list,array_sum_score_names)

#----- Save the Final_individual scores file for each "select_index_I"
save(FINAL_INDIVIDUAL_SCORES,file = paste0(results_dir,"/df_FINAL_INDIVIDUAL_SCORES_",exp_name,"_listIndex_",select_index_I,".RData"))
    

#-----------------------------------------------------------------------------
# PART 3: Regression model
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Run linear models: DEL and DUP

#----- 1. Regression model for DEL
str_covariates_lm_model_DEL <- paste0("Count_genes_within_list_DEL + Count_genes_not_in_list_DEL") # Add LM Covariates here
eval(parse(text = paste0("model_DEL =lm(",flag_phenotype_C," ~ ",str_covariates_lm_model_DEL,", data = FINAL_INDIVIDUAL_SCORES, na.action=na.omit)")))

# ----- save Model _DEL as RData file
save(model_DEL,file = paste0(results_dir,"/model_DEL_",exp_name,"_listIndex_",select_index_I,".RData"))

#----- 2. Regression model for DUP
str_covariates_lm_model_DUP <- paste0("Count_genes_within_list_DUP + Count_genes_not_in_list_DUP") # Add LM Covariates here
eval(parse(text = paste0("model_DUP =lm(",flag_phenotype_C," ~ ",str_covariates_lm_model_DUP,", data = FINAL_INDIVIDUAL_SCORES, na.action=na.omit)")))

# ----- save Model _DUP as RData file
save(model_DUP,file = paste0(results_dir,"/model_DUP_",exp_name,"_listIndex_",select_index_I,".RData"))



#-----------------------------------------------------------------------------
# END
#-----------------------------------------------------------------------------
#

