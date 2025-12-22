###################################################################
#### Generate a compound bargraph for DE protein counts per K-cluster 
#### for selected comparison groups
###################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Define directories
dir0 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dir1 <- file.path(dir0, "Kclust_all")
dirDATA <- file.path(dir1, "significant_csv")
output_dir <- file.path(dir1, "kCluster_bargraph")
if (!dir.exists(output_dir)) dir.create(output_dir)

# List only relevant CSV files
csv_files <- list.files(dirDATA, pattern = "^GFP_vs_.*\\.csv$", full.names = TRUE)

# Initialize list to store cluster counts
all_counts <- list()

# Process each CSV file
for (csv_file in csv_files) {
  df <- read.csv(csv_file)
  
  group_name <- basename(csv_file) %>%
    sub("_significant_peptides\\.csv$", "", .) %>%
    sub("_significant_proteins\\.csv$", "", .)
  
  type <- ifelse(grepl("proteins", csv_file, ignore.case = TRUE), "proteins", "peptides")
  
  if (nrow(df) < 1 || !"Cluster" %in% colnames(df)) {
    cluster_counts <- data.frame(
      Cluster = "NoData",
      Count = 1,
      Group = group_name,
      Type = type
    )
  } else {
    cluster_counts <- as.data.frame(table(df$Cluster))
    colnames(cluster_counts) <- c("Cluster", "Count")
    cluster_counts$Group <- group_name
    cluster_counts$Type <- type
  }
  
  all_counts[[length(all_counts) + 1]] <- cluster_counts
}

# Combine all into one data frame
combined_counts <- bind_rows(all_counts)

# Reassign clusters to K1–K7
unique_clusters <- sort(as.numeric(as.character(unique(combined_counts$Cluster[combined_counts$Cluster != "NoData"]))))

cluster_map <- setNames(paste0("K", seq_along(unique_clusters)), unique_clusters)
combined_counts$Cluster <- recode(combined_counts$Cluster, !!!cluster_map)

# Add NoData if present
if ("NoData" %in% combined_counts$Cluster) {
  combined_counts$Cluster[combined_counts$Cluster == "NoData"] <- "NoData"
}

# Set factor levels for Cluster
cluster_levels <- c(paste0("K", seq_along(unique_clusters)), "NoData")
combined_counts$Cluster <- factor(combined_counts$Cluster, levels = cluster_levels)

# Create full grid to ensure all combinations are present
full_grid <- expand.grid(
  Cluster = cluster_levels,
  Group = unique(combined_counts$Group),
  Type = unique(combined_counts$Type),
  stringsAsFactors = FALSE
)

complete_counts <- full_grid %>%
  left_join(combined_counts, by = c("Cluster", "Group", "Type")) %>%
  mutate(Count = ifelse(is.na(Count), 0, Count),
         FillGroup = ifelse(Cluster == "NoData", "NoData", Group))

# Define darker color palette
group_levels <- unique(complete_counts$FillGroup)
n_groups <- length(setdiff(group_levels, "NoData"))
dark_colors <- RColorBrewer::brewer.pal(n = max(3, min(8, n_groups)), name = "Set2")
group_colors <- setNames(c(dark_colors[1:n_groups], "gray30"), c(setdiff(group_levels, "NoData"), "NoData"))

# Plot function
# Remove NoData from plotting data and color mapping
plot_cluster_counts <- function(data, type_label, output_file) {
  plot_data <- filter(data, Type == type_label & Cluster != "NoData")
  
  # Recalculate group colors without NoData
  plot_groups <- unique(plot_data$FillGroup)
  plot_colors <- group_colors[plot_groups]
  
  p <- ggplot(plot_data, aes(x = Cluster, y = Count, fill = FillGroup)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = plot_colors) +
    theme_minimal() +
    labs(title = paste("Count - DE", type_label),
         x = "Cluster", y = "Count", fill = "Group")
  
  ggsave(file.path(output_dir, output_file), plot = p, width = 10, height = 6)
}

# Generate plots
plot_cluster_counts(complete_counts, "proteins", "Proteins_kCluster_Barplot.pdf")
plot_cluster_counts(complete_counts, "peptides", "Peptides_kCluster_Barplot.pdf")

# Export final data frame to CSV
write.csv(complete_counts, file = file.path(output_dir, "Final_kCluster_Counts.csv"), row.names = FALSE)
