#################################################################
#### This script is to generate K-means cluster based clustering
#### After differential expression analysis (saved SEs)
#################################################################


################### Step 1: Setup

library(ComplexHeatmap)
library(RColorBrewer)
library(matrixStats)
library(SummarizedExperiment)
library(tidyverse)
library(magick)

# Directories
dir0 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dirK <- file.path(dir0, "Kclust_all")
if (!dir.exists(dirK)) dir.create(dirK)

dirDE <- file.path(dir0, "DE_results")

# Load SE objects
load(file.path(dirDE, "de_results_bpca_list.RData"))
load(file.path(dirDE, "de_results_man_list.RData"))

all_se_lists <- list(
  bpca = get("de_results_bpca_list"),
  man = get("de_results_man_list")
)

################### Step 2: Heatmap Parameters
normalization_type <- "row_avg_centered" # Options: "row_avg_centered", "row_median_centered", "z_score_standardized"
scale_before_kmeans = FALSE # Choose FALSE as the matrix is log2 and centered
use_kmeans <- TRUE
k_clusters <- 7
k_seed = 997 # previously 123 
col_limit <- 4
clustering_distance <- "pearson" # Options: "euclidean", "pearson", "spearman", "manhattan"
clustering_method <- "complete" # Options: "complete", "ward.D2", "average"

################### Step 3: Heatmap Function 

custom_cluster_heatmap <- function(se, name, out_dir,
                                   normalization_type,
                                   use_kmeans,
                                   k_clusters,
                                   col_limit,
                                   clustering_distance,
                                   clustering_method,
                                   row_km_repeats = 100,
                                   export_matrix = TRUE,
                                   scale_before_kmeans = FALSE # [NEW]
) {
  expr <- assay(se)
  row_data <- rowData(se)
  col_data <- colData(se)
  
  if (!"significant" %in% names(row_data)) {
    message(" Skipping ", name, ": 'significant' column not found.")
    return(NULL)
  }
  
  sig_rows <- which(row_data$significant == TRUE)
  if (length(sig_rows) == 0) {
    message(" Skipping ", name, ": no significant features.")
    return(NULL)
  }
  
  mat <- expr[sig_rows, , drop = FALSE]
  
  # Normalization
  if (normalization_type == "row_avg_centered") {
    mat <- sweep(mat, 1, rowMeans(mat, na.rm = TRUE))
  } else if (normalization_type == "row_median_centered") {
    mat <- sweep(mat, 1, matrixStats::rowMedians(mat, na.rm = TRUE))
  } else if (normalization_type == "z_score_standardized") {
    mat <- sweep(mat, 1, rowMeans(mat, na.rm = TRUE))
    mat <- sweep(mat, 1, matrixStats::rowSds(mat, na.rm = TRUE), "/")
  }
  
  # Export normalized matrix
  if (export_matrix) {
    write.csv(mat, file.path(out_dir, paste0("normalized_matrix_", name, ".csv")))
  }
  
  # Optional: Scale data before K-means
  # Scaling before K-means should be avoided as the values are already log2 scale and centered
  if (use_kmeans && nrow(mat) >= k_clusters &&
      scale_before_kmeans &&
      !(normalization_type == "z_score_standardized") &&
      clustering_distance %in% c("euclidean", "spearman")) {
    message(" Scaling data for K-means: ", name)
    mat <- t(scale(t(mat)))
  }
  
  # K-means clustering
  cluster_annotation <- NULL
  row_ha <- NULL
  if (use_kmeans && nrow(mat) >= k_clusters) {
    message(" Running robust K-means for ", name)
    best_km <- NULL
    best_tot_withinss <- Inf
    
    for (i in 1:row_km_repeats) {
      set.seed(997 + i)
      km_try <- tryCatch(
        kmeans(mat, centers = k_clusters, iter.max = 100),
        error = function(e) NULL
      )
      if (!is.null(km_try) && km_try$tot.withinss < best_tot_withinss) {
        best_km <- km_try
        best_tot_withinss <- km_try$tot.withinss
      }
    }
    
    if (!is.null(best_km)) {
      km <- best_km
      mat <- mat[order(km$cluster), ]
      cluster_annotation <- km$cluster[order(km$cluster)]
      
      cluster_df <- data.frame(
        Feature = rownames(mat),
        Cluster = cluster_annotation
      )
      write.csv(cluster_df, file.path(out_dir, paste0("clusters_", name, ".csv")), row.names = FALSE)
      
      write.csv(data.frame(
        TotalWithinSS = km$tot.withinss,
        Betweenss = km$betweenss,
        Iterations = km$iter
      ), file.path(out_dir, paste0("kmeans_metrics_", name, ".csv")), row.names = FALSE)
      
      cluster_colors <- RColorBrewer::brewer.pal(max(3, length(unique(cluster_annotation))), "Set3")
      names(cluster_colors) <- as.character(sort(unique(cluster_annotation)))
      
      row_ha <- rowAnnotation(
        Cluster = factor(cluster_annotation),
        col = list(Cluster = cluster_colors),
        show_annotation_name = TRUE,
        annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),  # Label "Cluster"
        gp = gpar(fontsize = 14)  # Cluster number font size? doesn't look it works
      )
      
      
    } else {
      message(" K-means failed for ", name, ": no valid clustering result.")
    }
  } else {
    message(" Skipping K-means for ", name, ": insufficient rows or disabled.")
  }
  
  # Column annotation
  sample_groups <- col_data$condition
  group_colors <- RColorBrewer::brewer.pal(length(unique(sample_groups)), "Set1")
  names(group_colors) <- unique(sample_groups)
  
  col_ha <- HeatmapAnnotation(
    Condition = sample_groups,
    col = list(Condition = group_colors),
    annotation_legend_param = list(
      Condition = list(title = "Condition",
                       at = unique(sample_groups),
                       labels = unique(sample_groups))
    )
  )
  
  # Title and legend label
  title_text <- paste("K-means based clustering (", normalization_type, ")", sep = "")
  value_label <- if (normalization_type == "z_score_standardized") "Z-score" else "Log2 LFQ"
  
  # Heatmap
  heatmap_obj <- Heatmap(
    mat,
    name = value_label,
    column_title = title_text,
    column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "darkblue"),
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 14),  # Increase column label font size
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_columns = clustering_distance,
    clustering_method_columns = clustering_method,
    clustering_distance_rows = clustering_distance,
    clustering_method_rows = clustering_method,
    row_split = if (!is.null(cluster_annotation)) factor(cluster_annotation) else NULL,
    top_annotation = col_ha,
    left_annotation = row_ha,
    heatmap_legend_param = list(
      legend_height = unit(4, "cm"),
      title_gp = gpar(fontsize = 12, fontface = "bold"),  # Legend title font
      labels_gp = gpar(fontsize = 12)  # Legend labels font ? doesn't look it works
    ),
    col = circlize::colorRamp2(c(-col_limit, 0, col_limit), c("midnightblue", "white", "firebrick"))
  )
  
  # Save to PDF
  pdf(file.path(out_dir, paste0("heatmap_", name, "_", normalization_type, ".pdf")),
      width = 8, height = 10)
  draw(heatmap_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}


################### Step 4: Run Heatmap Generation 

for (list_name in names(all_se_lists)) {
  current_list <- all_se_lists[[list_name]]
  
  for (name in names(current_list)) {
    message("Processing: ", list_name, " - ", name)
    se <- current_list[[name]]
    
    custom_cluster_heatmap(
      se = se,
      name = name,
      out_dir = dirK,
      normalization_type = normalization_type,
      use_kmeans = use_kmeans,
      k_clusters = k_clusters,
      col_limit = col_limit,
      clustering_distance = clustering_distance,
      clustering_method = clustering_method,
      row_km_repeats = 100, # [NEW]
      export_matrix = TRUE, # [NEW]
      scale_before_kmeans = scale_before_kmeans
    )
  }
}
