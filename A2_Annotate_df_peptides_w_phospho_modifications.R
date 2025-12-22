library(tidyverse)

###### The Directories

dir250527 = "/media/akira/LeGris4h/250527_FragPipeAnalystR_1.0.5_MaxLFQ_Proteome_Prashant_Reorganized/20250527"

mod_peptide_file = "combined_modified_peptide.tsv"

dir00 <- "/media/akira/LeGris4h/250609_FragPipeAnalystR_1.1.0_MaxLFQ_250527"

# Load peptide data
peptide_df <- read_tsv(file.path(dir250527, mod_peptide_file))
colnames(peptide_df) <- gsub(" ", "_", colnames(peptide_df))

#### Update the obsolete gene symbols TTC26 and C5orf51
peptide_df$Gene[peptide_df$Gene == "TTC26"] = "IFT56"
peptide_df$Gene[peptide_df$Gene == "C5orf51"] = "RIMOC1"

#### Make a small subset df that only contains key ID columns without quantitative value columns
peptide_subset <- peptide_df %>%
  select(Peptide_Sequence, Modified_Sequence, Gene, Protein_ID, Protein_Description)

######## Annotate the df_bpca_peptides.csv -- this df includes LFQ values in replicates and is after bpca missing value imputation
### Select key columns and logical columns for DE
df_peptide_csv = "df_bpca_peptides.csv"
df_pep = read.csv(file.path(dir00, df_peptide_csv))

#### Select pY containing sequences in "peptide_tagged"
### Peptides with phosphorylated residues
pY_peptides <- peptide_subset %>%
  filter(str_detect(Modified_Sequence, "Y\\[79\\.9663\\]"))

pT_peptides <- peptide_subset %>%
  filter(str_detect(Modified_Sequence, "T\\[79\\.9663\\]"))

pS_peptides <- peptide_subset %>%
  filter(str_detect(Modified_Sequence, "S\\[79\\.9663\\]"))

## Step 1: Summarize each phospho dataset
# IMPORTANT: Handle multiple entries by summarizing with any() to ensure that if any instance of a peptide is phosphorylated on a residue, it is marked TRUE
pY_summary <- pY_peptides %>%
  group_by(Peptide_Sequence) %>%
  summarise(pTyr = any(TRUE), .groups = 'drop')

pT_summary <- pT_peptides %>%
  group_by(Peptide_Sequence) %>%
  summarise(pThr = any(TRUE), .groups = 'drop')

pS_summary <- pS_peptides %>%
  group_by(Peptide_Sequence) %>%
  summarise(pSer = any(TRUE), .groups = 'drop')

## Step 2: Start with your main dataframe, df_pep
annotated_df <- df_pep %>%
  left_join(pY_summary, by = "Peptide_Sequence") %>%
  left_join(pT_summary, by = "Peptide_Sequence") %>%
  left_join(pS_summary, by = "Peptide_Sequence")

## Step 3: Replace NAs with FALSE
annotated_df <- annotated_df %>%
  mutate(across(c(pTyr, pThr, pSer), ~replace_na(.x, FALSE)))

## Step 4: Create a combination column
# The Combination column encodes the presence (1) or absence (0) of each phosphorylation type in the order: pTyr, pThr, pSer
annotated_df <- annotated_df %>%
  mutate(Combination = paste0(as.integer(pTyr), as.integer(pThr), as.integer(pSer)))

## Step 5: Take average of LFQ Intensity between the replicates
# Identify columns ending with "_MaxLFQ_Intensity"
intensity_cols <- grep("_MaxLFQ_Intensity$", names(annotated_df), value = TRUE)

# Extract base group names (assuming replicates are named like Group1_1_MaxLFQ_Intensity)
group_names <- unique(gsub("_[12]_MaxLFQ_Intensity$", "", intensity_cols))

# Loop through each group and compute the average of the two replicates
for (group in group_names) {
  rep1 <- paste0(group, "_1_MaxLFQ_Intensity")
  rep2 <- paste0(group, "_2_MaxLFQ_Intensity")
  avg_col <- paste0(group, "_Avg_MaxLFQ_Intensity")
  
  # Compute row-wise mean and add to dataframe
  annotated_df[[avg_col]] <- rowMeans(annotated_df[, c(rep1, rep2)], na.rm = TRUE)
}

## Save the modified df_bpca_peptide dataframe
write.csv(annotated_df, file.path(dir00, "consensus/df_bpca_peptide_phosphoMAP.csv"), row.names = F)



########## Box Plots and ANOVA

# Function to prepare log2-transformed long-format data
prepare_log2_data <- function(df, filter_condition = NULL) {
  if (!is.null(filter_condition)) {
    df <- df %>% filter(!!rlang::parse_expr(filter_condition))
    if (nrow(df) == 0) return(NULL)
  }
  avg_cols <- grep("_Avg_MaxLFQ_Intensity$", names(df), value = TRUE)
  df[avg_cols] <- lapply(df[avg_cols], function(x) ifelse(x == 0, NA, x))
  log2_df <- log2(df[avg_cols])
  log2_long <- pivot_longer(log2_df, cols = everything(), names_to = "Group", values_to = "Log2_Intensity")
  log2_long$Group <- gsub("_Avg_MaxLFQ_Intensity", "", log2_long$Group)
  return(na.omit(log2_long))
}

# Function to plot boxplot
plot_log2_boxplot <- function(log2_long, output_file, plot_title) {
  if (is.null(log2_long)) {
    message("No data to plot.")
    return()
  }
  pdf(output_file)
  print(
    ggplot(log2_long, aes(x = Group, y = Log2_Intensity, fill = Group)) +
      geom_boxplot(color = "black") +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() +
      labs(title = plot_title,
           x = "Experimental Group",
           y = "Log2 Avg_MaxLFQ_Intensity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  )
  dev.off()
}

# Function to run ANOVA and Tukey HSD
run_anova_tukey <- function(log2_long, output_csv) {
  if (is.null(log2_long)) {
    message("No data for statistical analysis.")
    return()
  }
  anova_result <- aov(Log2_Intensity ~ Group, data = log2_long)
  print(summary(anova_result))
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
  write.csv(as.data.frame(tukey_result$Group), file = output_csv)
}

### Use the functions defined above
# Set output directory
output_dir <- file.path(dir00, "consensus")

# All peptides
log2_all <- prepare_log2_data(annotated_df)
plot_log2_boxplot(log2_all,
                  file.path(output_dir, "boxplot_avg_maxLFQ_all.pdf"),
                  "Box Plot: Log2 Avg MaxLFQ Intensity for all peptides")
run_anova_tukey(log2_all,
                file.path(output_dir, "anova_tukey_all.csv"))

# Filtered peptides (Combination == '100')
log2_pY <- prepare_log2_data(annotated_df, "Combination == '100'")
plot_log2_boxplot(log2_pY,
                  file.path(output_dir, "boxplot_avg_maxLFQ_pY_notpT_notpS.pdf"),
                  "Box Plot: Log2 Avg MaxLFQ Intensity for the Phospho Pattern pY-notpT-notpS")
run_anova_tukey(log2_pY,
                file.path(output_dir, "anova_tukey_pY.csv"))

