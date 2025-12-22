###########
#  Pairwise Fisher's Exact Tests with FDR Comparing Clusters within each Group
###########

# Load required packages
# install.packages("dplyr")
library(dplyr)

# --- Load and filter data ---
# Note: Using your local file path
dirDATA <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527/Kclust_all/kCluster_bargraph"
data <- "Final_kCluster_Counts.csv"
df <- read.csv(file.path(dirDATA, data))
df_filtered <- df[df$Cluster != "NoData", ]

# --- Create output directories ---
dir.create(file.path(dirDATA, "peptide_fisher_pairwise_per_group"), showWarnings = FALSE)
dir.create(file.path(dirDATA, "protein_fisher_pairwise_per_group"), showWarnings = FALSE)

# --- Define groups and types to loop through ---
unique_groups <- unique(df_filtered$Group)

# --- Run Analysis for each Group and Type ---
for (group_name in unique_groups) {
  for (type_name in unique_types) {
    
    subset_data <- subset(df_filtered, Group == group_name & Type == type_name)
    
    if (nrow(subset_data) >= 2) {
      cat(paste("Processing:", type_name, "in Group", group_name, "...\n"))
      
      total_count_for_subset <- sum(subset_data$Count)
      cluster_pairs <- combn(subset_data$Cluster, 2, simplify = FALSE)
      results_for_subset <- list()
      
      for (pair in cluster_pairs) {
        cluster1_name <- pair[1]
        cluster2_name <- pair[2]
        
        count1 <- subset_data$Count[subset_data$Cluster == cluster1_name]
        count2 <- subset_data$Count[subset_data$Cluster == cluster2_name]
        
        contingency_table <- matrix(c(
          count1, total_count_for_subset - count1,
          count2, total_count_for_subset - count2
        ), nrow = 2, byrow = TRUE)
        
        test_result <- fisher.test(contingency_table)
        
        results_for_subset[[length(results_for_subset) + 1]] <- data.frame(
          Cluster1 = cluster1_name,
          Cluster2 = cluster2_name,
          p_value = test_result$p.value
        )
      }
      
      if (length(results_for_subset) > 0) {
        final_results_df <- do.call(rbind, results_for_subset)
        final_results_df$p_adjusted_fdr <- p.adjust(final_results_df$p_value, method = "fdr")
        
        output_dir <- ifelse(type_name == "peptides", 
                             file.path(dirDATA, "peptide_fisher_pairwise_per_group"), 
                             file.path(dirDATA, "protein_fisher_pairwise_per_group"))
        
        output_file <- file.path(output_dir, paste0("Group_", group_name, "_pairwise_fisher.csv"))
        write.csv(final_results_df, output_file, row.names = FALSE)
      }
    }
  }
}

cat("\nAnalysis 2 (Fisher's Test) complete!\n")