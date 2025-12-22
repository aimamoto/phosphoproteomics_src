###########
#  Pairwise Fisher's Exact Tests with FDR Comparing Groups within each Cluster
###########

# Load required packages
library(dplyr)

# --- Load and filter data ---
dirDATA <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527/Kclust_all/kCluster_bargraph"
data <- "Final_kCluster_Counts.csv"
df <- read.csv(file.path(dirDATA, data))
df_filtered <- df[df$Cluster != "NoData", ]

# --- Create output directories ---
dir.create(file.path(dirDATA, "peptide_fisher_pairwise_per_cluster"), showWarnings = FALSE)
dir.create(file.path(dirDATA, "protein_fisher_pairwise_per_cluster"), showWarnings = FALSE)

# --- Define clusters and types to loop through ---
unique_clusters <- unique(df_filtered$Cluster)
unique_types <- c("peptides", "proteins")

# --- Run Analysis for each Cluster and Type ---
for (clust in unique_clusters) {
  for (type_name in unique_types) {
    
    subset_data <- subset(df_filtered, Cluster == clust & Type == type_name)
    
    if (nrow(subset_data) >= 2) {
      cat(paste("Processing:", type_name, "in Cluster", clust, "...\n"))
      
      # Calculate the total count for the current cluster/type
      total_count_for_subset <- sum(subset_data$Count)
      
      group_pairs <- combn(subset_data$Group, 2, simplify = FALSE)
      results_for_subset <- list()
      
      for (pair in group_pairs) {
        group1_name <- pair[1]
        group2_name <- pair[2]
        
        count1 <- subset_data$Count[subset_data$Group == group1_name]
        count2 <- subset_data$Count[subset_data$Group == group2_name]
        
        # Create the 2x2 contingency table
        contingency_table <- matrix(c(
          count1, total_count_for_subset - count1,
          count2, total_count_for_subset - count2
        ), nrow = 2, byrow = TRUE)
        
        # Run Fisher's test
        test_result <- fisher.test(contingency_table)
        
        results_for_subset[[length(results_for_subset) + 1]] <- data.frame(
          Group1 = group1_name,
          Group2 = group2_name,
          p_value = test_result$p.value
        )
      }
      
      if (length(results_for_subset) > 0) {
        final_results_df <- do.call(rbind, results_for_subset)
        final_results_df$p_adjusted_fdr <- p.adjust(final_results_df$p_value, method = "fdr")
        
        output_dir <- ifelse(type_name == "peptides", 
                             file.path(dirDATA, "peptide_fisher_pairwise_per_cluster"), 
                             file.path(dirDATA, "protein_fisher_pairwise_per_cluster"))
        
        output_file <- file.path(output_dir, paste0("Cluster_", clust, "_pairwise_fisher.csv"))
        write.csv(final_results_df, output_file, row.names = FALSE)
      }
    }
  }
}

cat("\nAnalysis 1 (Fisher's Test) complete!\n")