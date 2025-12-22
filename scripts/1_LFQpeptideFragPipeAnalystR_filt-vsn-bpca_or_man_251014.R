############ Proteome Analysis 
## Prashant has revised the FragPipe peptide and protein files for 2 bioreplicates on 250527
## As FragPipe-Analyst web or shiny app versions have serious issues, I am using FragPipeAnalystR, a dedicated R package (https://github.com/Nesvilab/FragPipeAnalystR)

#library(DEP) ## Add additional functions. See https://github.com/arnesmits/DEP NOTE: Attaching the package DEP globally causes errors in some places, for example, with 'bpcn' in 'pcaMethods'. Call it locally only for specific functions by "DEP::"
library(FragPipeAnalystR)
library(SummarizedExperiment)
library(tidyverse)

##### The Directories and data files
dir0 = "/media/akira/argentee/proteome/250527_FragPipeAnalystR_1.0.5_MaxLFQ_Proteome_Prashant_Reorganized/20250527"

pep = file.path(dir0,"combined_peptide.tsv")
pro = file.path(dir0,"combined_protein.tsv")
expANNO = file.path(dir0, "experiment_annotation_manually_edited3.tsv")

dir1 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
if(!dir.exists(dir1)) dir.create(dir1)

dirSE = file.path(dir1, "saved_SEs")
dirPCA = file.path(dir1, "PCA")
dirDE = file.path(dir1, "DE_results")
dirVOL = file.path(dir1, "volcano")

if(!dir.exists(dirSE)) dir.create(dirSE)
if(!dir.exists(dirPCA)) dir.create(dirPCA)

if(!dir.exists(dirDE)) dir.create(dirDE)
if(!dir.exists(dirVOL)) dir.create(dirVOL)

##### Generate Summarized Objects (SEs) using the function of FragPipeAnalystR
se_pep = make_se_from_files(pep, expANNO, type = "LFQ", level = "peptide")
se_pro = make_se_from_files(pro, expANNO, type = "LFQ", level = "protein")

## It is a minor issue, but two gene names TTC26 and C5orf51 are obsolete and need to be updated to IFT56 and RIMOC1, respectively
rowData(se_pep)$Gene[rowData(se_pep)$Gene == "TTC26"] <- "IFT56"
rowData(se_pro)$Gene[rowData(se_pro)$Gene == "TTC26"] <- "IFT56"

rowData(se_pep)$Gene[rowData(se_pep)$Gene == "C5orf51"] <- "RIMOC1"
rowData(se_pro)$Gene[rowData(se_pro)$Gene == "C5orf51"] <- "RIMOC1"

##################### Missing Value Imputation (MVI)
## Filter missing values: first -- but encountered error in 'left_join()'
#filt_se_pep = DEP::filter_missval(se_pep, thr = 1)
#filt_se_pro = DEP::filter_missval(se_pro, thr = 1)

## Custom prefilter for Missing Value Reduction
# This filter will remove the rows the have a proportion of NAs in each row greater than a value of threshold
# NOTE: The lower the threshold value is, the more stringent the prefilter will become.
# CAUTION: Sample groups that have more MVs/NAs will be affected more than groups with fewer NAs. Control group often have greater numbers of NAs.

filterSE <- function(se, threshold, assay_name = 1) {
  # 1. Access the specified assay matrix from the object
  assay_matrix <- assay(se, assay_name)
  
  # 2. Calculate the proportion of NAs for each row (feature)
  na_proportion <- rowSums(is.na(assay_matrix)) / ncol(assay_matrix)
  
  # 3. Create a logical vector of which rows to KEEP
  rows_to_keep <- na_proportion <= threshold
  
  # 4. Subset the original SummarizedExperiment object
  filtered_se <- se[rows_to_keep, ]
  
  # 5. Return the new, smaller SummarizedExperiment object
  return(filtered_se)
}

# For moderate stringency, try a threshold of 0.5 (50%)
filt50_se_pep = filterSE(se_pep, threshold = 0.5)
filt50_se_pro = filterSE(se_pro, threshold = 0.5)

## Optional: Variance Stabilizing Normalization -- Calling the package DEP locally. See the note at the beginning of the script
vsn_50_se_pep = DEP::normalize_vsn(filt50_se_pep)
vsn_50_se_pro = DEP::normalize_vsn(filt50_se_pro)

#library(hexbin)
#DEP::meanSdPlot(se_pep)
#DEP::meanSdPlot(vsn_50_se_pep)

## MVI options: c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "man", "min", "zero", "mixed", "nbavg")
## Choose Bayesian PCA missing value estimation "bpca", But it takes more time than the option below
imputed_bpca_se_pep = impute(se_pep, fun = "bpca")
imputed_bpca_se_pro = impute(se_pro, fun = "bpca")

imputed_bpca_50_pep = impute(filt50_se_pep, fun = "bpca")
imputed_bpca_50_pro = impute(filt50_se_pro, fun = "bpca")

imputed_bpca_vsn_50_pep = impute(vsn_50_se_pep, fun = "bpca")
imputed_bpca_vsn_50_pro = impute(vsn_50_se_pro, fun = "bpca")

## OPTION: Perseus Type MVI: Select "man" Two parameters are set by default as "scale = 0.3, shift = 1.8"
imputed_man_se_pep = impute(se_pep, fun = "man")
imputed_man_se_pro = impute(se_pro, fun = "man")

imputed_man_50_pep = impute(filt50_se_pep, fun = "man")
imputed_man_50_pro = impute(filt50_se_pro, fun = "man")

imputed_man_vsn_50_pep = impute(vsn_50_se_pep, fun = "man")
imputed_man_vsn_50_pro = impute(vsn_50_se_pro, fun = "man")


#### Define a function to convert rowData to tibble and clean column names
process_row_data <- function(se_object) {
  df <- as_tibble(rowData(se_object))
  colnames(df) <- gsub("\\.", "_", colnames(df))
  return(df)
}

### 1. Define the SE objects in a named list for bpca MVI groups
bpca_se_list <- list(
  bpca_result_pep = imputed_bpca_se_pep,
  bpca_result_pro = imputed_bpca_se_pro,
  bpca_50_result_pep = imputed_bpca_50_pep,
  bpca_50_result_pro = imputed_bpca_50_pro,
  bpca_vsn_50_result_pep = imputed_bpca_vsn_50_pep,
  bpca_vsn_50_result_pro = imputed_bpca_vsn_50_pro
)

# Apply the function to each item in the list
df_bpca_list <- lapply(bpca_se_list, process_row_data)

### 2. Define the SE objects in a named list for Perseus-type MVI groups
man_se_list <- list(
  man_result_pep = imputed_man_se_pep,
  man_result_pro = imputed_man_se_pro,
  man_50_result_pep = imputed_man_50_pep,
  man_50_result_pro = imputed_man_50_pro,
  man_vsn_50_result_pep = imputed_man_vsn_50_pep,
  man_vsn_50_result_pro = imputed_man_vsn_50_pro
)

# Apply the function to each item in the list
df_man_list <- lapply(man_se_list, process_row_data)

### 3. Select specific columns for each dataframe in the list

# Function for peptide-level data frames
select_pep_columns <- function(df) {
  df %>%
    select(
      "Peptide_Sequence",
      "Protein_ID",
      "Gene",
      "Description",
      "imputed",
      ends_with("_Intensity")
    )
}

# Function for protein-level data frames
select_pro_columns <- function(df) {
  df %>%
    select(
      "Protein_ID",
      "Gene",
      "Description",
      "imputed",
      ends_with("_Intensity")
    )
}

## 3.1 for bpca list
# Apply the appropriate function based on the name
df_bpca_list <- lapply(names(df_bpca_list), function(name) {
  df <- df_bpca_list[[name]]
  if (grepl("_pep$", name)) {
    select_pep_columns(df)
  } else {
    select_pro_columns(df)
  }
})

# Preserve names
names(df_bpca_list) <- names(bpca_se_list)

## 3.2 for man list
df_man_list <- lapply(names(df_man_list), function(name) {
  df <- df_man_list[[name]]
  if (grepl("_pep$", name)) {
    select_pep_columns(df)
  } else {
    select_pro_columns(df)
  }
})

# Preserve names
names(df_man_list) <- names(man_se_list)

## 3.3 Save each data frame with a custom filename
mapply(function(df, name) {
  # Customize the filename (e.g., remove "df_" prefix if needed)
  filename <- paste0("df_", name, ".csv")
  filepath <- file.path(dirSE, filename)
  write.csv(df, filepath, row.names = FALSE)
}, df_bpca_list, names(df_bpca_list))

mapply(function(df, name) {
  # Customize the filename (e.g., remove "df_" prefix if needed)
  filename <- paste0("df_", name, ".csv")
  filepath <- file.path(dirSE, filename)
  write.csv(df, filepath, row.names = FALSE)
}, df_man_list, names(df_man_list))

################# PCA Plots

# Loop through the SE list and save PCA plots
mapply(function(se_obj, name) {
  # Customize the filename
  filename <- paste0("PCA_", name, ".pdf")
  filepath <- file.path(dirPCA, filename)
  
  # Save the plot
  pdf(filepath)
  print(plot_pca(se_obj))
  dev.off()
}, bpca_se_list, names(bpca_se_list))

mapply(function(se_obj, name) {
  # Customize the filename
  filename <- paste0("PCA_", name, ".pdf")
  filepath <- file.path(dirPCA, filename)
  
  # Save the plot
  pdf(filepath)
  print(plot_pca(se_obj))
  dev.off()
}, man_se_list, names(man_se_list))

################# Testing Differential Expression 

#### For bpca MVI SEs
# Run 'test_limma' and 'add_rejections' for each SE object
de_results_bpca_list <- lapply(bpca_se_list, function(se_obj) {
  de_result <- test_limma(se_obj, type = "all")
  de_result <- add_rejections(de_result)
  return(de_result)
})

# Name the list for easy access
names(de_results_bpca_list) <- names(bpca_se_list)

#### For Perseus MVI SEs
# Run test_limma and add_rejections for each SE object
de_results_man_list <- lapply(man_se_list, function(se_obj) {
  de_result <- test_limma(se_obj, type = "all")
  de_result <- add_rejections(de_result)
  return(de_result)
})

# Name the list for easy access
names(de_results_man_list) <- names(man_se_list)

#### Save the differential expression results

save(de_results_bpca_list, file = file.path(dirDE, "de_results_bpca_list.RData"))
save(de_results_man_list, file = file.path(dirDE, "de_results_man_list.RData"))

#### Let's save the SEs used for DE analysis as well
save(bpca_se_list, file = file.path(dirSE, "bpca_se_list.RData"))
save(man_se_list, file = file.path(dirSE, "man_se_list.RData"))

################# Volcano Plots
#### Selected sets of 1:1 comparisons
contrast_set1 = c("GFP_vs_P8E2", "GFP_vs_WT", "P8E2_vs_WT")
contrast_set2 = c("P8E2_vs_P8E2_T338I", "P8E2_vs_P8E2_T338M", "P8E2_T338I_vs_P8E2_T338M")
contrast_set3 = c("WT_vs_WT_T338I", "WT_vs_WT_T338M", "WT_T338I_vs_WT_T338M")

all_contrasts <- list(contrast_set1, contrast_set2, contrast_set3)

### 1. For 'de_results_bpca_list'
for (name in names(de_results_bpca_list)) {
  de_result <- de_results_bpca_list[[name]]
  
  for (contrast_set in all_contrasts) {
    for (contrast in contrast_set) {
      # Define output file name using SE name
      filename <- paste0("volcano_", name, "_", contrast, ".pdf")
      filepath <- file.path(dirVOL, filename)
      
      # Determine if it's protein-level (ends with "_pro")
      is_protein <- grepl("_pro$", name)
      
      # Save the volcano plot
      pdf(filepath)
      print(
        if (is_protein) {
          plot_volcano(de_result, contrast, name_col = "Gene") +
            ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
            ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40")
        } else {
          plot_volcano(de_result, contrast) +
            ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
            ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40")
        }
      )
      dev.off()
    }
  }
}

### 2. For 'de_results_man_list'
for (name in names(de_results_man_list)) {
  de_result <- de_results_man_list[[name]]
  
  for (contrast_set in all_contrasts) {
    for (contrast in contrast_set) {
      # Define output file name using SE name
      filename <- paste0("volcano_", name, "_", contrast, ".pdf")
      filepath <- file.path(dirVOL, filename)
      
      # Determine if it's protein-level (ends with "_pro")
      is_protein <- grepl("_pro$", name)
      
      # Save the volcano plot
      pdf(filepath)
      print(
        if (is_protein) {
          plot_volcano(de_result, contrast, name_col = "Gene") +
            ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
            ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40")
        } else {
          plot_volcano(de_result, contrast) +
            ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
            ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40")
        }
      )
      dev.off()
    }
  }
}


################# Extract rowData as dataframe from the SE Objects tested for DE
# Function to extract and clean rowData
extract_row_data <- function(se_obj) {
  df <- as_tibble(rowData(se_obj))
  colnames(df) <- gsub("\\.", "_", colnames(df))
  return(df)
}

### 1. for the 'de_result_bpca_list'
# Apply to each SE object in the list
df_DE_bpca_list <- lapply(de_results_bpca_list, extract_row_data)

# Name the list for easy access
names(df_DE_bpca_list) <- names(de_results_bpca_list)

### 2. for the 'de_result_man_list'
# Apply to each SE object in the list
df_DE_man_list <- lapply(de_results_man_list, extract_row_data)

# Name the list for easy access
names(df_DE_man_list) <- names(de_results_man_list)

## Select the limma results and a few key columns
# Function for peptide-level data frames
simplify_pep_df <- function(df) {
  df %>%
    select(
      "Peptide_Sequence",
      "Protein_ID",
      "Gene",
      "Description",
      "imputed",
      "significant",
      ends_with("_p_val"),
      ends_with("_p_adj"),
      ends_with("_diff"),
      ends_with("_significant")
    )
}

# Function for protein-level data frames
simplify_pro_df <- function(df) {
  df %>%
    select(
      "Protein_ID",
      "Gene",
      "Description",
      "imputed",
      "significant",
      ends_with("_p_val"),
      ends_with("_p_adj"),
      ends_with("_diff"),
      ends_with("_significant")
    )
}

### 1. For df_DE_bpca_list
# Apply the appropriate function based on the name
df_DE_bpca_simple_list <- lapply(names(df_DE_bpca_list), function(name) {
  df <- df_DE_bpca_list[[name]]
  if (grepl("_pep$", name)) {
    simplify_pep_df(df)
  } else {
    simplify_pro_df(df)
  }
})

# Preserve names
names(df_DE_bpca_simple_list) <- names(de_results_bpca_list)

### 2. For df_DE_man_list
# Apply the appropriate function based on the name
df_DE_man_simple_list <- lapply(names(df_DE_man_list), function(name) {
  df <- df_DE_man_list[[name]]
  if (grepl("_pep$", name)) {
    simplify_pep_df(df)
  } else {
    simplify_pro_df(df)
  }
})

# Preserve names
names(df_DE_man_simple_list) <- names(de_results_man_list)


################# Count the number of "significant" peptides or proteins

# Function to summarize significant and total counts
summarize_significance <- function(df) {
  c(
    Significant_Count = sum(df$significant, na.rm = TRUE),
    Total_Count = sum(!is.na(df$significant))
  )
}

### 1. For 'df_DE_bpca_simple_list'
# Apply to each data frame in the list
summary_bpca_list <- lapply(df_DE_bpca_simple_list, summarize_significance)

# Convert to a data frame for easier viewing
summary_bpca <- do.call(rbind, summary_bpca_list)
summary_bpca <- as.data.frame(summary_bpca)
summary_bpca$SE_Name <- rownames(summary_bpca)
rownames(summary_bpca) <- NULL
summary_bpca <- summary_bpca[, c("SE_Name", "Significant_Count", "Total_Count")]

### 2. For 'df_DE_man_simple_list'
# Apply to each data frame in the list
summary_man_list <- lapply(df_DE_man_simple_list, summarize_significance)

# Convert to a data frame for easier viewing
summary_man <- do.call(rbind, summary_man_list)
summary_man <- as.data.frame(summary_man)
summary_man$SE_Name <- rownames(summary_man)
rownames(summary_man) <- NULL
summary_man <- summary_man[, c("SE_Name", "Significant_Count", "Total_Count")]

################ Save the dataframe csv files from list
write.csv(summary_bpca, file.path(dirDE, "summary_counts_bpca_FragPipeAnalystR.csv"), row.names = F)
write.csv(summary_man, file.path(dirDE, "summary_counts_Perseus_FragPipeAnalystR.csv"), row.names = F)

### Save each data frame in df_DE_list
# 1. bpca
mapply(function(df, name) {
  filename <- paste0("DE_", name, ".csv")
  filepath <- file.path(dirDE, filename)
  write.csv(df, filepath, row.names = FALSE)
}, df_DE_bpca_simple_list, names(df_DE_bpca_simple_list))

# 2. Perseus
mapply(function(df, name) {
  filename <- paste0("DE_", name, ".csv")
  filepath <- file.path(dirDE, filename)
  write.csv(df, filepath, row.names = FALSE)
}, df_DE_man_simple_list, names(df_DE_man_simple_list))

