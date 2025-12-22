###############  
library(dplyr)

##### 1. Identify the proteins flagged as DE in comparisons 
### DE is judged as TRUE or FALSE in the output file
### Change the data directory dir_data for your local environment!

dir0 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
data_dir = file.path(dir0, "DE_results")
filename_pep = "DE_bpca_vsn_50_result_pep.csv"
filename_pro = "DE_bpca_vsn_50_result_pro.csv"

clust_dir = file.path(dir0, "Kclust_all")
filename_clust_pep = "clusters_bpca_vsn_50_result_pep.csv"
filename_clust_pro = "clusters_bpca_vsn_50_result_pro.csv"

dirOUT = file.path(clust_dir, "significant_csv")
if (!dir.exists(dirOUT)) {
  dir.create(dirOUT)
}

### Read data for DE peptide and proteins and filter only for significant
dfpep = read.csv(file.path(data_dir, filename_pep))
dfpro = read.csv(file.path(data_dir, filename_pro))

dfpep = dfpep %>%
  filter(significant == TRUE)
dfpro = dfpro %>%
  filter(significant == TRUE)

### Read K-means clusters
dfclst_pep = read.csv(file.path(clust_dir, filename_clust_pep))
dfclst_pro = read.csv(file.path(clust_dir, filename_clust_pro))

# Remove the proteinID before and including the underscore, as the Feature column has a string combined proteinID and sequence like "P36578_AAAAAAALQAKSDEK"
dfclst_pep$Feature = gsub("^[^_]*_", "", dfclst_pep$Feature)

# Add a new column Cluster in dfpep
dfpep = dfpep %>%
  left_join(dfclst_pep, by = c("Peptide_Sequence" = "Feature"))

# Add a new column Cluster in dfpro
dfpro = dfpro %>%
  left_join(dfclst_pro, by = c("Protein_ID" = "Feature"))

### These are the name of all columns with "_significant" that indicate logical calls as TRUE or FALSE
compGRPs_pro = names(dfpro[grep("_significant", names(dfpro))])
compGRPs_pep = names(dfpep[grep("_significant", names(dfpep))])

prot = c("Protein_ID", "Gene", "Cluster")
pept = c("Peptide_Sequence", "Protein_ID", "Gene", "Cluster")


#### Define a new function to identify significant proteins 
findsigPROs = function(dfpro, compGRPs_pro) { 
  sigPRO_list = list()
  # Create a logical vector for significance 
  for (grp in compGRPs_pro) {
  sigPROs = dfpro %>%
    filter(dfpro[[grp]] == TRUE) %>%
    select(all_of(prot))
  sigPRO_list[[grp]] = sigPROs
  }
  return(sigPRO_list)
}

#### Define a new function to identify significant peptides 
findsigPEPs = function(dfpep, compGRPs_pep) { 
  sigPEP_list = list()
  # Create a logical vector for significance 
  for (grp in compGRPs_pep) {
    sigPEPs = dfpep %>%
      filter(dfpep[[grp]] == TRUE) %>%
      select(all_of(pept))
    sigPEP_list[[grp]] = sigPEPs
  }
  return(sigPEP_list)
}

# Identify significant proteins and peptides for each comparison group
sigPRO_list = findsigPROs(dfpro, compGRPs_pro)
sigPEP_list = findsigPEPs(dfpep, compGRPs_pep)

# Save the list of significant proteins for each comparison group
for (grp in compGRPs_pro) {
  write.csv(sigPRO_list[[grp]],file.path(dirOUT, paste0(grp,"_proteins.csv")), row.names = F)
}

for (grp in compGRPs_pep) {
  write.csv(sigPEP_list[[grp]],file.path(dirOUT, paste0(grp,"_peptides.csv")), row.names = F)
}
