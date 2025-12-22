###########
#  Pairwise Chi-Square Tests Comparing Clusters within each Group
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
dir.create(file.path(dirDATA, "peptide_chisq_pairwise_per_group"), showWarnings = FALSE)
dir.create(file.path(dirDATA, "protein_chisq_pairwise_per_group"), showWarnings = FALSE)

# --- Define groups and types to loop through ---
unique_groups <- unique(df_filtered$Group)
unique_types <- c("peptides", "proteins")

# --- Run Analysis for each Group and Type ---
for (group_name in unique_groups) {
  for (type_name in unique_types) {
    
    # Subset the data for the current Group and Type
    subset_data <- subset(df_filtered, Group == group_name & Type == type_name)
    
    # Proceed only if there are at least two clusters to compare
    if (nrow(subset_data) >= 2) {
      
      cat(paste("Processing:", type_name, "in Group", group_name, "...\n"))
      
      # Step 1: Create all unique 2-cluster combinations
      cluster_pairs <- combn(subset_data$Cluster, 2, simplify = FALSE)
      
      # Initialize a list to store results for this subset
      results_for_subset <- list()
      
      # Step 2: Loop through each pair and perform a Chi-squared test
      for (pair in cluster_pairs) {
        cluster1_name <- pair[1]
        cluster2_name <- pair[2]
        
        # Get the counts for the two clusters
        counts <- c(
          subset_data$Count[subset_data$Cluster == cluster1_name],
          subset_data$Count[subset_data$Cluster == cluster2_name]
        )
        
        # Run the test (only if both counts are not zero)
        if (sum(counts) > 0) {
          test_result <- chisq.test(counts)
          
          # Store the result
          results_for_subset[[length(results_for_subset) + 1]] <- data.frame(
            Cluster1 = cluster1_name,
            Cluster2 = cluster2_name,
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
                             file.path(dirDATA, "peptide_chisq_pairwise_per_group"), 
                             file.path(dirDATA, "protein_chisq_pairwise_per_group"))
        
        output_file <- file.path(output_dir, paste0("Group_", group_name, "_pairwise.csv"))
        write.csv(final_results_df, output_file, row.names = FALSE)
      }
    }
  }
}

cat("\nAll analyses complete!\n")