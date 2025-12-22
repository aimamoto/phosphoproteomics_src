#################
# Let's make a DE result with K-means clustering numbers
#################

dir0 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"

#### Read DE_result files for both peptide and protein
dirDE = file.path(dir0, "DE_results")
DE_list = c(
  "DE_bpca_vsn_50_result_pep.csv",
  "DE_bpca_vsn_50_result_pro.csv"
  )

dfDE_list = lapply(file.path(dirDE, DE_list), read.csv)

# Remove columns ending with "_p_val" from each dataframe in the list
dfDE_list <- lapply(dfDE_list, function(df) {
  df[ , !grepl("_p_val$", names(df))]
})

# Keep the rows that are TRUE for the column "significant"
dfDE_list_sig = lapply(dfDE_list, function(df) df[df$significant == TRUE,])

#### Read K-Cluster files for both peptide and protein
dirKclust = file.path(dir0, "Kclust_all")
Kclust_list = c(
  "clusters_bpca_vsn_50_result_pep.csv", 
  "clusters_bpca_vsn_50_result_pro.csv"
  )

dfKclust_list = lapply(file.path(dirKclust, Kclust_list), read.csv)

# Remove the string before "_" in the Feature column of "clusters_bpca_vsn_50_result_pep.csv"
dfKclust_list[[1]]$Feature = gsub(".*?_", "", dfKclust_list[[1]]$Feature)

# For dfDE_list_sig[[1]]: match on Peptide_Sequence
dfDE_list_sig[[1]] <- merge(
  dfDE_list_sig[[1]],
  dfKclust_list[[1]][, c("Feature", "Cluster")],
  by.x = "Peptide_Sequence",
  by.y = "Feature",
  all.x = TRUE
)

# For dfDE_list_sig[[2]]: match on Protein_ID
dfDE_list_sig[[2]] <- merge(
  dfDE_list_sig[[2]],
  dfKclust_list[[2]][, c("Feature", "Cluster")],
  by.x = "Protein_ID",
  by.y = "Feature",
  all.x = TRUE
)

names(dfDE_list_sig) = c(
  "DE_bpca_vsn_50_Kclust_pep", 
  "DE_bpca_vsn_50_Kclust_pro"
  )

##### Separate DE pep or pro by Cluster
# Create an empty list to store the new dataframes
cluster_split_list <- list()

# Loop over each dataframe in dfDE_list_sig
for (name in names(dfDE_list_sig)) {
  df <- dfDE_list_sig[[name]]
  
  # Determine prefix based on name
  prefix <- if (grepl("pep", name)) "pepK" else "proK"
  
  # Split by Cluster and assign names
  for (k in 1:7) {
    df_k <- df[df$Cluster == k, ]
    if (nrow(df_k) > 0) {
      cluster_split_list[[paste0(prefix, k)]] <- df_k
    }
  }
}

dirDE_Kclust = file.path(dirKclust, "DE_by_Kclusters")
if (!dir.exists(dirDE_Kclust)) dir.create(dirDE_Kclust)

# Save each dataframe in dfDE_list_sig to a CSV file
Map(function(df, name) {
  write.csv(df, file = file.path(dirDE_Kclust, paste0(name, ".csv")), row.names = FALSE)
}, dfDE_list_sig, names(dfDE_list_sig))

# Save each dataframe in cluster_split_list to a CSV file
Map(function(df, name) {
  write.csv(df, file = file.path(dirDE_Kclust, paste0(name, ".csv")), row.names = FALSE)
}, cluster_split_list, names(cluster_split_list))
