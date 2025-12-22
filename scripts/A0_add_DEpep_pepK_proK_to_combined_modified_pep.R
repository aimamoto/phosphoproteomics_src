#########################################################################
#### Modify the combined_modified_peptide.tsv
#### Add K-means clusters (protein cluster and peptide cluster) 
#### Add DE for peptides and proteins
#### Keep Peptide_Sequence, Modified_Sequence, Protein_ID, Gene, 
#########################################################################

#### Step 1: Simplify combined_modified_peptide.tsv
dirD = "/media/akira/argentee/proteome/250527_FragPipeAnalystR_1.0.5_MaxLFQ_Proteome_Prashant_Reorganized/20250527"

library(tidyverse)

# Read main df
dfpepMOD = read_tsv(file.path(dirD, "combined_modified_peptide.tsv"))
# Replace empty space with underscore
colnames(dfpepMOD) = gsub(" ", "_", colnames(dfpepMOD))

## It is a minor issue, but two gene names TTC26 and C5orf51 are obsolete and need to be updated to IFT56 and RIMOC1, respectively
dfpepMOD$Gene[dfpepMOD$Gene == "TTC26"] <- "IFT56"
dfpepMOD$Gene[dfpepMOD$Gene == "C5orf51"] <- "RIMOC1"

dfpepMOD_simple = dfpepMOD %>%
  select(
    Peptide_Sequence, 
    Modified_Sequence, 
    Protein_ID,
    Gene,
    ends_with("_MaxLFQ_Intensity")
         )

# Read DE results by peptides and proteins both with K clusters
dir0 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dirK = file.path(dir0, "Kclust_all")
dirDE = file.path(dirK, "DE_by_Kclusters")
pepKDE = read.csv(file.path(dirDE, "DE_bpca_vsn_50_Kclust_pep.csv"))
proKDE = read.csv(file.path(dirDE, "DE_bpca_vsn_50_Kclust_pro.csv"))

# Add new columns to dfpepMOD_simple for Cluster by the Peptide_Sequence
# Select only the relevant columns from pepKDE
pepKDE_subset <- pepKDE %>%
  select(Peptide_Sequence, contains("_significant"), Cluster)

# Join the dataframes on Peptide_Sequence
dfpepMOD_simple <- dfpepMOD_simple %>%
  left_join(pepKDE_subset, by = "Peptide_Sequence") %>%
  drop_na() 

dfpepMOD_simple = dfpepMOD_simple %>%
  mutate(Cluster = paste0("pepK", Cluster))

dfpepMOD_simple <- dfpepMOD_simple %>%
  rename(pep_Cluster = Cluster)

# Add a new column to dfpepMOD_simple for Cluster 
dfpepMOD_simple <- dfpepMOD_simple %>%
  left_join(proKDE %>% select(Protein_ID, Cluster), by = "Protein_ID") %>%
  rename(pro_Cluster = Cluster)

# Change the cluster number to proK#, but leave NA as NA
dfpepMOD_simple <- dfpepMOD_simple %>%
  mutate(pro_Cluster = if_else(
    is.na(pro_Cluster),
    NA_character_,
    paste0("proK", pro_Cluster)
  ))

# Change order by pro_Cluster
dfpepMOD_simple <- dfpepMOD_simple %>%
  arrange(pro_Cluster)

write.csv(dfpepMOD_simple, file.path(dirDE, "combined_modified_DE_Kclust.csv"), row.names = F)

# Optional, but just select the rows for pro_Cluster == proK4 and sort for GFP_vs_P8E2
dfpepMOD_simple_proK4 = dfpepMOD_simple %>%
  filter(pro_Cluster == "proK4") %>%
  arrange(desc(GFP_vs_P8E2_significant))

write.csv(dfpepMOD_simple_proK4, file.path(dirDE, "combined_modified_DE_proK4_only.csv"), row.names = F)

