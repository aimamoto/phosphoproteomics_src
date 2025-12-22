# =======================
# CompareCluster Enrichment Analysis Script
# =======================

# ---- Load Libraries ----
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(enrichplot)

########### The Directories and Input Data ----
dir0 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"

dirDE = file.path(dir0, "DE_results")
dirK <- file.path(dir0, "Kclust_all")

clustDA <- "clusters_bpca_vsn_50_result_pro.csv"
deDA = "DE_bpca_vsn_50_result_pro.csv"

dfclust <- read.csv(file.path(dirK, clustDA))

dfDE = read.csv(file.path(dirDE, deDA))
dfDEsig = dfDE %>%
  filter(significant == TRUE)

dirOUT <- file.path(dirK, "enrich_by_kClusters_per_DE_GFPvP8orWT")
if (!dir.exists(dirOUT)) dir.create(dirOUT)

########### Map UniProt IDs to Entrez and Symbol
dfclust_mapped <- dfclust %>%
  mutate(
    ENTREZID = mapIds(org.Hs.eg.db,
                      keys = Feature,
                      column = "ENTREZID",
                      keytype = "UNIPROT",
                      multiVals = "first"),
    SYMBOL = mapIds(org.Hs.eg.db,
                    keys = Feature,
                    column = "SYMBOL",
                    keytype = "UNIPROT",
                    multiVals = "first")
  )

#### To focus on GFPvP8E2 and GFPvWT, dfclust_mapped need to be attached with 'GFP_vs_P8E2_significant' and 'GFP_vs_WT_sigificant' (and other '_significant' columns) from dfDEsig by dfDEsig$Protein_ID vs dfclust_mapped$Feature

  # Step 1: Select only the columns with "_significant" in their names, plus the key column
  signif_cols <- grep("_significant", names(dfDEsig), value = TRUE)
  dfDEsig_subset <- dfDEsig[, c("Protein_ID", signif_cols)]
    
  # Step 2: Merge with dfclust_mapped using Feature and Protein_ID as keys
  dfclust_merged <- merge(dfclust_mapped, dfDEsig_subset, by.x = "Feature", by.y = "Protein_ID", all.x = TRUE)

#### Filter dfclust_merged by 'GFP_vs_P8E2_significant' and 'GFP_vs_WT_sigificant'

dfclust_GFPvP8 = dfclust_merged %>%
  filter(GFP_vs_P8E2_significant == TRUE)
dfclust_GFPvWT = dfclust_merged %>%
  filter(GFP_vs_WT_significant == TRUE)

# ---- Prepare Gene List for compareCluster ----
gene_list_P8 <- dfclust_GFPvP8 %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(Cluster) %>%
  summarise(genes = list(unique(ENTREZID))) %>%
  deframe()

gene_list_WT <- dfclust_GFPvWT %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(Cluster) %>%
  summarise(genes = list(unique(ENTREZID))) %>%
  deframe()

##### ---- Filter clusters that have less than 5 proteins (arbitrary cutoff)
gene_list_P8 <- gene_list_P8[sapply(gene_list_P8, length) >= 10]
gene_list_WT <- gene_list_WT[sapply(gene_list_WT, length) >= 10]


##### ---- Run compareCluster for Multiple Ontologies ----

# GO-BP
compare_go_bp_P8 <- compareCluster(geneCluster = gene_list_P8,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pvalueCutoff = 0.05)

if (!is.null(compare_go_bp_P8)) {
  simple_go_bp_P8 <- simplify(compare_go_bp_P8)
} else {
  message("No enrichment found for GO:MF in P8 clusters.")
}

compare_go_bp_WT <- compareCluster(geneCluster = gene_list_WT,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "BP",
                                   pvalueCutoff = 0.05)
simple_go_bp_WT = simplify(compare_go_bp_WT)

# GO-MF
compare_go_mf_P8 <- compareCluster(geneCluster = gene_list_P8,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "MF",
                                   pvalueCutoff = 0.05)


if (!is.null(compare_go_mf_P8)) {
  simple_go_mf_P8 <- simplify(compare_go_mf_P8)
} else {
  message("No enrichment found for GO:MF in P8 clusters.")
}

compare_go_mf_WT <- compareCluster(geneCluster = gene_list_WT,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "MF",
                                   pvalueCutoff = 0.05)
simple_go_mf_WT = simplify(compare_go_mf_WT)

# GO-CC
compare_go_cc_P8 <- compareCluster(geneCluster = gene_list_P8,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "CC",
                                   pvalueCutoff = 0.05)

if (!is.null(compare_go_cc_P8)) {
  simple_go_cc_P8 <- simplify(compare_go_cc_P8)
} else {
  message("No enrichment found for GO:MF in P8 clusters.")
}

compare_go_cc_WT <- compareCluster(geneCluster = gene_list_WT,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "CC",
                                   pvalueCutoff = 0.05)
simple_go_cc_WT = simplify(compare_go_cc_WT)

# KEGG
compare_kegg_P8 <- compareCluster(geneCluster = gene_list_P8,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.05)
compare_kegg_WT <- compareCluster(geneCluster = gene_list_WT,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.05)

# MSigDB Hallmark
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

compare_hallmark_P8 <- compareCluster(geneCluster = gene_list_P8,
                                   fun = "enricher",
                                   TERM2GENE = msig_hallmark %>% select(gs_name, ncbi_gene),
                                   pvalueCutoff = 0.05)

compare_hallmark_WT <- compareCluster(geneCluster = gene_list_WT,
                                   fun = "enricher",
                                   TERM2GENE = msig_hallmark %>% select(gs_name, ncbi_gene),
                                   pvalueCutoff = 0.05)

# WikiPathways
msig_wiki <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")

compare_wiki_P8 <- compareCluster(geneCluster = gene_list_P8,
                               fun = "enricher",
                               TERM2GENE = msig_wiki %>% select(gs_name, ncbi_gene),
                               pvalueCutoff = 0.05)

compare_wiki_WT <- compareCluster(geneCluster = gene_list_WT,
                                  fun = "enricher",
                                  TERM2GENE = msig_wiki %>% select(gs_name, ncbi_gene),
                                  pvalueCutoff = 0.05)


############# ---- Save Dotplots ----

###### simplified GO
# GO_BP
pdf(file.path(dirOUT, "compareCluster_simpleGO_BP_P8.pdf"), width = 6, height = 10) # Cluster dotplot is on a smaller canvas size as fewer terms appear in this group
print(
  dotplot(simple_go_bp_P8, showCategory = 8) + 
    ggtitle("Simplified GO-BP by cluster for GFP vs P8E2") +
    theme(axis.text.y = element_text(size = 12), 
          plot.title = element_text(size = 12, face = "bold") # Plot title
    )
  )
dev.off()

pdf(file.path(dirOUT, "compareCluster_simpleGO_BP_WT.pdf"), width = 10, height = 15)
print(
  dotplot(simple_go_bp_WT, showCategory = 8) +
    ggtitle("Simplified GO-BP by cluster for GFP vs WT") +
    theme(axis.text.y = element_text(size = 15), 
          
          # GO term labels
          axis.text.x = element_text(size = 14), # Cluster labels
          axis.title.x = element_text(size = 16), # "Cluster" title
          axis.title.y = element_text(size = 14), # "Legends" title
          plot.title = element_text(size = 18, face = "bold"), # Plot title
          legend.title = element_text(size = 13), # Legend title
          legend.text = element_text(size = 11) # Legend items
          )
)
dev.off()

# GO_MF
if (!is.null(compare_go_mf_P8)) {
  pdf(file.path(dirOUT, "compareCluster_simpleGO_MF_P8.pdf"), width = 10, height = 15)
  print(dotplot(simple_go_mf_P8, showCategory = 8) + ggtitle("Simplified GO Molecular Function by cluster for GFP vs P8E2"))
  dev.off()
} else {
  message("Skipping GO:MF P8 plot — no enrichment results.")
}

pdf(file.path(dirOUT, "compareCluster_simpleGO_MF_WT.pdf"), width = 10, height = 15)
print(dotplot(simple_go_mf_WT, showCategory = 8) + ggtitle("Simplified GO Molecular Function by cluster for GFP vs WT"))
dev.off()

# GO_CC
pdf(file.path(dirOUT, "compareCluster_simpleGO_CC_P8.pdf"), width = 10, height = 15)
print(dotplot(simple_go_cc_P8, showCategory = 8) + ggtitle("Simplified GO Cellular Component by cluster for GFP vs P8E2"))
dev.off()

pdf(file.path(dirOUT, "compareCluster_simpleGO_CC_WT.pdf"), width = 10, height = 15)
print(dotplot(simple_go_cc_WT, showCategory = 8) + ggtitle("Simplified GO Cellular Component by cluster for GFP vs WT"))
dev.off()

# KEGG
if (!is.null(compare_kegg_P8)) {
  pdf(file.path(dirOUT, "compareCluster_KEGG_P8.pdf"), width = 10, height = 15)
  print(dotplot(compare_kegg_P8, showCategory = 8) + ggtitle("KEGG Enrichment by cluster for GFP vs P8E2"))
  dev.off()
} else {
  message("Skipping KEGG P8 plot — no enrichment results.")
}

pdf(file.path(dirOUT, "compareCluster_KEGG_WT.pdf"), width = 10, height = 15)
print(dotplot(compare_kegg_WT, showCategory = 8) + ggtitle("KEGG Enrichment by cluster for GFP vs WT"))
dev.off()

# Hallmark
pdf(file.path(dirOUT, "compareCluster_Hallmark_P8.pdf"), width = 10, height = 15)
print(dotplot(compare_hallmark_P8, showCategory = 8) + ggtitle("Hallmark Enrichment by cluster for GFP vs P8E2"))
dev.off()

pdf(file.path(dirOUT, "compareCluster_Hallmark_WT.pdf"), width = 10, height = 15)
print(dotplot(compare_hallmark_WT, showCategory = 8) + ggtitle("Hallmark Enrichment by cluster for GFP vs WT"))
dev.off()

# WikiPathway
pdf(file.path(dirOUT, "compareCluster_WikiPathway_P8.pdf"), width = 10, height = 15)
print(dotplot(compare_wiki_P8, showCategory = 8) + ggtitle("WikiPathway Enrichment by cluster for GFP vs P8E2"))
dev.off()

pdf(file.path(dirOUT, "compareCluster_WikiPathway_WT.pdf"), width = 10, height = 15)
print(dotplot(compare_wiki_WT, showCategory = 8) + ggtitle("WikiPathway Enrichment by cluster for GFP vs WT"))
dev.off()


################## ---- Save Tables ----

# KEGG
write.csv(as.data.frame(compare_kegg_P8), file.path(dirOUT, "compareCluster_KEGG_P8.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_kegg_WT), file.path(dirOUT, "compareCluster_KEGG_WT.csv"), row.names = FALSE)

# Hallmark
write.csv(as.data.frame(compare_hallmark_P8), file.path(dirOUT, "compareCluster_Hallmark_P8.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_hallmark_WT), file.path(dirOUT, "compareCluster_Hallmark_WT.csv"), row.names = FALSE)

# WikiPathways
write.csv(as.data.frame(compare_wiki_P8), file.path(dirOUT, "compareCluster_WikiPathways_P8.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_wiki_WT), file.path(dirOUT, "compareCluster_WikiPathways_WT.csv"), row.names = FALSE)

# Simplified GO
write.csv(as.data.frame(simple_go_bp_P8), file.path(dirOUT, "compareCluster_simpleGO_BP_P8.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_bp_WT), file.path(dirOUT, "compareCluster_simpleGO_BP_WT.csv"), row.names = FALSE)

write.csv(as.data.frame(simple_go_mf_P8), file.path(dirOUT, "compareCluster_simpleGO_MF_P8.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_mf_WT), file.path(dirOUT, "compareCluster_simpleGO_MF_WT.csv"), row.names = FALSE)

write.csv(as.data.frame(simple_go_cc_P8), file.path(dirOUT, "compareCluster_simpleGO_CC_P8.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_cc_WT), file.path(dirOUT, "compareCluster_simpleGO_CC_WT.csv"), row.names = FALSE)
