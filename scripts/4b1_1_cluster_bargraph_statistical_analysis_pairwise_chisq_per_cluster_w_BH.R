###########
#  Manual Pairwise Chi-square Tests with FDR for each K-means cluster
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
dir.create(file.path(dirDATA, "peptide_chisq_pairwise_per_cluster"), showWarnings = FALSE)
dir.create(file.path(dirDATA, "protein_chisq_pairwise_per_cluster"), showWarnings = FALSE)

# --- Define clusters and types to loop through ---
unique_clusters <- unique(df_filtered$Cluster)
unique_types <- c("peptides", "proteins")

# --- Run Analysis for each Cluster and Type ---
for (clust in unique_clusters) {
  for (type_name in unique_types) {
    
    # Subset the data
    subset_data <- subset(df_filtered, Cluster == clust & Type == type_name)
    
    # Proceed only if there are at least two groups to compare
    if (nrow(subset_data) >= 2) {
      
      cat(paste("Processing:", type_name, "in Cluster", clust, "...\n"))
      
      # Step 1: Create all unique 2-group combinations
      group_pairs <- combn(subset_data$Group, 2, simplify = FALSE)
      
      # Initialize a list to store results for this subset
      results_for_subset <- list()
      
      # Step 2: Loop through each pair and perform a Chi-squared test
      for (pair in group_pairs) {
        group1_name <- pair[1]
        group2_name <- pair[2]
        
        # Get the counts for the two groups
        counts <- c(
          subset_data$Count[subset_data$Group == group1_name],
          subset_data$Count[subset_data$Group == group2_name]
        )
        
        # Run the test (only if both counts are not zero)
        if (sum(counts) > 0) {
          test_result <- chisq.test(counts)
          
          # Store the result
          results_for_subset[[length(results_for_subset) + 1]] <- data.frame(
            Group1 = group1_name,
            Group2 = group2_name,
            p_value = test_result$p.value
          )
        }
      }
      
      # Step 3: Combine results and adjust p-values
      if (length(results_for_subset) > 0) {
        final_results_df <- do.call(rbind, results_for_subset)
        final_results_df$p_adjusted_fdr <- p.adjust(final_results_df$p_value, method = "fdr")
        
        # Define the output path and save the file
        output_dir <- ifelse(type_name == "peptides", 
                             file.path(dirDATA, "peptide_chisq_pairwise_per_cluster"), 
                             file.path(dirDATA, "protein_chisq_pairwise_per_cluster"))
        
        output_file <- file.path(output_dir, paste0("Cluster_", clust, "_pairwise.csv"))
        write.csv(final_results_df, output_file, row.names = FALSE)
      }
    }
  }
}

cat("\nAll analyses complete!\n")