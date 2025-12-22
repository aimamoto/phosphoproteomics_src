### This script convert geneID (ENTREZID) to gene symbol

## NOTE: GO enrichments list genes with Symbols (will skip to PART II in this script)

############## The Directories
dir0 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"

dirdata = file.path(dir0, "Kclust_all")
dirdataP8WT = file.path(dirdata, "enrich_by_kClusters_per_DE_GFPvP8orWT")
dirdata338I338M = file.path(dirdata, "enrich_by_kClusters_per_DE_GFPvT338IorT338M")
dirdataP8338I338M = file.path(dirdata, "enrich_by_kClusters_per_DE_GFPvP8T338IorP8T338M")

###1: Define and create the new output directory ---
# Define the path for the new output directory
dirOUT = file.path(dirdata, "GO_geneSymbols_enrich_by_kClusters")
# Create the directory if it doesn't already exist
if (!dir.exists(dirOUT)) dir.create(dirOUT)
# -----------------------------------------------------------


library(org.Hs.eg.db)
library(AnnotationDbi)

################ PART I List of CSV files to process for geneID to symbol conversion
# Define the regular expression pattern to match the file names
# This pattern finds files containing "compareCluster_simpleGO_BP" OR "compareCluster_simpleGO_CC" and that end with ".csv"
file_pattern <- "compareCluster_simpleGO_(BP|CC|MF).*\\.csv$"

# List the matching files from both directories
# full.names = TRUE ensures you get the complete file path
files_P8WT <- list.files(path = dirdataP8WT, 
                         pattern = file_pattern, 
                         full.names = TRUE)

files_338I338M <- list.files(path = dirdata338I338M, 
                             pattern = file_pattern, 
                             full.names = TRUE)

files_P8338I338M <- list.files(path = dirdataP8338I338M, 
                               pattern = file_pattern, 
                               full.names = TRUE)

# Combine the two lists of file paths into a single vector
all_target_files <- c(files_P8WT, files_338I338M, files_P8338I338M)

# Custom function
extract_gene_info <- function(df, geneID, cluster_col, description_col) {
  results <- list()
  for (i in 1:nrow(df)) {
    gene_ids <- unlist(strsplit(df[i, geneID], "/"))
    gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    gene_names <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")
    temp_df <- data.frame(
      Cluster = df[i, cluster_col],
      Description = df[i, description_col],
      geneID = gene_ids,
      geneSymbol = gene_symbols,
      geneName = gene_names
    )
    results[[i]] <- temp_df
  }
  new_df <- do.call(rbind, results)
  return(new_df)
}

# Process each file
for (file in all_target_files) {
  df <- read.csv(file)
  new_df <- extract_gene_info(df, "geneID", "Cluster", "Description")
  
  # --- CHANGE 2: Create the output filename using the new directory ---
  # Get just the filename from the full path (e.g., "my_data.csv")
  input_filename <- basename(file)
  # Change the suffix of the base filename
  output_filename_base <- sub(".csv", "_by_Symbol.csv", input_filename)
  # Join the new directory path with the new filename
  output_file <- file.path(dirOUT, output_filename_base)
  # ------------------------------------------------------------------
  
  write.csv(new_df, output_file, row.names = FALSE)
}