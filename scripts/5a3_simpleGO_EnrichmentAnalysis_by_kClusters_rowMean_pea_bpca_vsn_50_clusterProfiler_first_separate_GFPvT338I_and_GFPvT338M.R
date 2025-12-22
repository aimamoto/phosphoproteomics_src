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

############ The Directories and Input Data
dir0 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"

dirDE = file.path(dir0, "DE_results")
dirK <- file.path(dir0, "Kclust_all")

clustDA <- "clusters_bpca_vsn_50_result_pro.csv"
deDA = "DE_bpca_vsn_50_result_pro.csv"

dfclust <- read.csv(file.path(dirK, clustDA))

dfDE = read.csv(file.path(dirDE, deDA))
dfDEsig = dfDE %>%
  filter(significant == TRUE)

# Map UniProt IDs to Entrez and Symbol
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

#### To focus on GFPvT338I and GFPvT338M, 'dfclust_mapped' need to be attached with 'GFP_vs_WT_T338I_significant' and 'GFP_vs_WT_T338M_sigificant' (and other '_significant' columns) from dfDEsig by dfDEsig$Protein_ID vs dfclust_mapped$Feature

  # Step 1: Select only the columns with "_significant" in their names, plus the key column
  signif_cols <- grep("_significant", names(dfDEsig), value = TRUE)
  dfDEsig_subset <- dfDEsig[, c("Protein_ID", signif_cols)]
    
  # Step 2: Merge with dfclust_mapped using Feature and Protein_ID as keys
  dfclust_merged <- merge(dfclust_mapped, dfDEsig_subset, by.x = "Feature", by.y = "Protein_ID", all.x = TRUE)

#### Filter dfclust_merged by 'GFP_vs_WT_T338I_significant' and 'GFP_vs_WT_T338M_sigificant'

dfclust_GFPvT338I = dfclust_merged %>%
  filter(GFP_vs_WT_T338I_significant == TRUE)
dfclust_GFPvT338M = dfclust_merged %>%
  filter(GFP_vs_WT_T338M_significant == TRUE)

# ---- Define Output Directory ----
dirOUT <- file.path(dirK, "enrich_by_kClusters_per_DE_GFPvT338IorT338M")
if (!dir.exists(dirOUT)) dir.create(dirOUT)

# ---- Prepare Gene List for compareCluster ----
gene_list_T338I <- dfclust_GFPvT338I %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(Cluster) %>%
  summarise(genes = list(unique(ENTREZID))) %>%
  deframe()

gene_list_T338M <- dfclust_GFPvT338M %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(Cluster) %>%
  summarise(genes = list(unique(ENTREZID))) %>%
  deframe()

##### ---- Filter clusters that have less than 5 proteins (arbitrary cutoff)
gene_list_T338I <- gene_list_T338I[sapply(gene_list_T338I, length) >= 10]
gene_list_T338M <- gene_list_T338M[sapply(gene_list_T338M, length) >= 10]


##### ---- Run compareCluster for Multiple Ontologies ----

# GO-BP
compare_go_bp_T338I <- compareCluster(geneCluster = gene_list_T338I,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pvalueCutoff = 0.05)

if (!is.null(compare_go_bp_T338I)) {
  simple_go_bp_T338I <- simplify(compare_go_bp_T338I)
} else {
  message("No enrichment found for GO:MF in T338I clusters.")
}

compare_go_bp_T338M <- compareCluster(geneCluster = gene_list_T338M,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "BP",
                                   pvalueCutoff = 0.05)
simple_go_bp_T338M = simplify(compare_go_bp_T338M)

# GO-MF
compare_go_mf_T338I <- compareCluster(geneCluster = gene_list_T338I,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "MF",
                                   pvalueCutoff = 0.05)


if (!is.null(compare_go_mf_T338I)) {
  simple_go_mf_T338I <- simplify(compare_go_mf_T338I)
} else {
  message("No enrichment found for GO:MF in T338I clusters.")
}

compare_go_mf_T338M <- compareCluster(geneCluster = gene_list_T338M,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "MF",
                                   pvalueCutoff = 0.05)
simple_go_mf_T338M = simplify(compare_go_mf_T338M)

# GO-CC
compare_go_cc_T338I <- compareCluster(geneCluster = gene_list_T338I,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "CC",
                                   pvalueCutoff = 0.05)

if (!is.null(compare_go_cc_T338I)) {
  simple_go_cc_T338I <- simplify(compare_go_cc_T338I)
} else {
  message("No enrichment found for GO:MF in T338I clusters.")
}

compare_go_cc_T338M <- compareCluster(geneCluster = gene_list_T338M,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   ont = "CC",
                                   pvalueCutoff = 0.05)
simple_go_cc_T338M = simplify(compare_go_cc_T338M)

# KEGG
compare_kegg_T338I <- compareCluster(geneCluster = gene_list_T338I,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.05)

compare_kegg_T338M <- compareCluster(geneCluster = gene_list_T338M,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.05)

# MSigDB Hallmark
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")

compare_hallmark_T338I <- compareCluster(geneCluster = gene_list_T338I,
                                   fun = "enricher",
                                   TERM2GENE = msig_hallmark %>% select(gs_name, ncbi_gene),
                                   pvalueCutoff = 0.05)

compare_hallmark_T338M <- compareCluster(geneCluster = gene_list_T338M,
                                   fun = "enricher",
                                   TERM2GENE = msig_hallmark %>% select(gs_name, ncbi_gene),
                                   pvalueCutoff = 0.05)

# WikiPathways
msig_wiki <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")

compare_wiki_T338I <- compareCluster(geneCluster = gene_list_T338I,
                               fun = "enricher",
                               TERM2GENE = msig_wiki %>% select(gs_name, ncbi_gene),
                               pvalueCutoff = 0.05)

compare_wiki_T338M <- compareCluster(geneCluster = gene_list_T338M,
                                  fun = "enricher",
                                  TERM2GENE = msig_wiki %>% select(gs_name, ncbi_gene),
                                  pvalueCutoff = 0.05)


############# ---- Save Dotplots ----

###### simplified GO
# GO_BP
pdf(file.path(dirOUT, "compareCluster_simpleGO_BP_T338I.pdf"), width = 10, height = 15)
print(
  dotplot(simple_go_bp_T338I, showCategory = 8) + 
    ggtitle("Simplified GO-BP by cluster for GFP vs WT_T338I") +
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

pdf(file.path(dirOUT, "compareCluster_simpleGO_BP_T338M.pdf"), width = 10, height = 15)
print(
  dotplot(simple_go_bp_T338M, showCategory = 8) +
    ggtitle("Simplified GO-BP by cluster for GFP vs WT_T338M") +
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
if (!is.null(compare_go_mf_T338I)) {
  pdf(file.path(dirOUT, "compareCluster_simpleGO_MF_T338I.pdf"), width = 10, height = 15)
  print(
    dotplot(simple_go_mf_T338I, showCategory = 8) + 
      ggtitle("Simplified GO-MF by cluster for GFP vs WT_T338I")) +
    theme(axis.text.y = element_text(size = 15), 
          
          # GO term labels
          axis.text.x = element_text(size = 14), # Cluster labels
          axis.title.x = element_text(size = 16), # "Cluster" title
          axis.title.y = element_text(size = 14), # "Legends" title
          plot.title = element_text(size = 18, face = "bold"), # Plot title
          legend.title = element_text(size = 13), # Legend title
          legend.text = element_text(size = 11) # Legend items
          )
  dev.off()
} else {
  message("Skipping GO:MF T338I plot — no enrichment results.")
}

pdf(file.path(dirOUT, "compareCluster_simpleGO_MF_T338M.pdf"), width = 10, height = 15)
print(dotplot(simple_go_mf_T338M, showCategory = 8) + 
        ggtitle("Simplified GO-MF by cluster for GFP vs WT_T338M")) +
  theme(axis.text.y = element_text(size = 15), 
        
        # GO term labels
        axis.text.x = element_text(size = 14), # Cluster labels
        axis.title.x = element_text(size = 16), # "Cluster" title
        axis.title.y = element_text(size = 14), # "Legends" title
        plot.title = element_text(size = 18, face = "bold"), # Plot title
        legend.title = element_text(size = 13), # Legend title
        legend.text = element_text(size = 11) # Legend items
  )
dev.off()

# GO_CC
pdf(file.path(dirOUT, "compareCluster_simpleGO_CC_T338I.pdf"), width = 10, height = 15)
print(dotplot(simple_go_cc_T338I, showCategory = 8) + 
        ggtitle("Simplified GO-CC by cluster for GFP vs WT_T338I"))  +
  theme(axis.text.y = element_text(size = 15), 
        
        # GO term labels
        axis.text.x = element_text(size = 14), # Cluster labels
        axis.title.x = element_text(size = 16), # "Cluster" title
        axis.title.y = element_text(size = 14), # "Legends" title
        plot.title = element_text(size = 18, face = "bold"), # Plot title
        legend.title = element_text(size = 13), # Legend title
        legend.text = element_text(size = 11) # Legend items
  )
dev.off()

pdf(file.path(dirOUT, "compareCluster_simpleGO_CC_T338M.pdf"), width = 10, height = 15)
print(dotplot(simple_go_cc_T338M, showCategory = 8) +
      ggtitle("Simplified GO-CC by cluster for GFP vs WT_T338M"))  +
  theme(axis.text.y = element_text(size = 15), 
        
        # GO term labels
        axis.text.x = element_text(size = 14), # Cluster labels
        axis.title.x = element_text(size = 16), # "Cluster" title
        axis.title.y = element_text(size = 14), # "Legends" title
        plot.title = element_text(size = 18, face = "bold"), # Plot title
        legend.title = element_text(size = 13), # Legend title
        legend.text = element_text(size = 11) # Legend items
  )
dev.off()

# KEGG
if (!is.null(compare_kegg_T338I)) {
  pdf(file.path(dirOUT, "compareCluster_KEGG_T338I.pdf"), width = 10, height = 15)
  print(dotplot(compare_kegg_T338I, showCategory = 8) +
        ggtitle("KEGG Enrichment by cluster for GFP vs WT_T338I"))  +
    theme(axis.text.y = element_text(size = 15), 
          
          # GO term labels
          axis.text.x = element_text(size = 14), # Cluster labels
          axis.title.x = element_text(size = 16), # "Cluster" title
          axis.title.y = element_text(size = 14), # "Legends" title
          plot.title = element_text(size = 18, face = "bold"), # Plot title
          legend.title = element_text(size = 13), # Legend title
          legend.text = element_text(size = 11) # Legend items
    )
  dev.off()
} else {
  message("Skipping KEGG T338I plot — no enrichment results.")
}

pdf(file.path(dirOUT, "compareCluster_KEGG_T338M.pdf"), width = 10, height = 15)
print(dotplot(compare_kegg_T338M, showCategory = 8) + 
        ggtitle("KEGG Enrichment by cluster for GFP vs WT_T338M"))  + 
  
  theme(axis.text.y = element_text(size = 15), 
        
        # GO term labels
        axis.text.x = element_text(size = 14), # Cluster labels
        axis.title.x = element_text(size = 16), # "Cluster" title
        axis.title.y = element_text(size = 14), # "Legends" title
        plot.title = element_text(size = 18, face = "bold"), # Plot title
        legend.title = element_text(size = 13), # Legend title
        legend.text = element_text(size = 11) # Legend items
  )
dev.off()

# Hallmark
pdf(file.path(dirOUT, "compareCluster_Hallmark_T338I.pdf"), width = 10, height = 15)
print(dotplot(compare_hallmark_T338I, showCategory = 8) + ggtitle("Hallmark Enrichment by cluster for GFP vs WT_T338I"))
dev.off()

pdf(file.path(dirOUT, "compareCluster_Hallmark_T338M.pdf"), width = 10, height = 15)
print(dotplot(compare_hallmark_T338M, showCategory = 8) + ggtitle("Hallmark Enrichment by cluster for GFP vs WT_T338M"))
dev.off()

# WikiPathway
pdf(file.path(dirOUT, "compareCluster_WikiPathway_T338I.pdf"), width = 10, height = 15)
print(dotplot(compare_wiki_T338I, showCategory = 8) + ggtitle("WikiPathway Enrichment by cluster for GFP vs WT_T338I"))
dev.off()

pdf(file.path(dirOUT, "compareCluster_WikiPathway_T338M.pdf"), width = 10, height = 15)
print(dotplot(compare_wiki_T338M, showCategory = 8) + ggtitle("WikiPathway Enrichment by cluster for GFP vs WT_T338M"))
dev.off()


################## ---- Save Tables ----

# KEGG
write.csv(as.data.frame(compare_kegg_T338I), file.path(dirOUT, "compareCluster_KEGG_T338I.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_kegg_T338M), file.path(dirOUT, "compareCluster_KEGG_T338M.csv"), row.names = FALSE)

# Hallmark
write.csv(as.data.frame(compare_hallmark_T338I), file.path(dirOUT, "compareCluster_Hallmark_T338I.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_hallmark_T338M), file.path(dirOUT, "compareCluster_Hallmark_T338M.csv"), row.names = FALSE)

# WikiPathways
write.csv(as.data.frame(compare_wiki_T338I), file.path(dirOUT, "compareCluster_WikiPathways_T338I.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_wiki_T338M), file.path(dirOUT, "compareCluster_WikiPathways_T338M.csv"), row.names = FALSE)

# Simplified GO
write.csv(as.data.frame(simple_go_bp_T338I), file.path(dirOUT, "compareCluster_simpleGO_BP_T338I.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_bp_T338M), file.path(dirOUT, "compareCluster_simpleGO_BP_T338M.csv"), row.names = FALSE)

write.csv(as.data.frame(simple_go_mf_T338I), file.path(dirOUT, "compareCluster_simpleGO_MF_T338I.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_mf_T338M), file.path(dirOUT, "compareCluster_simpleGO_MF_T338M.csv"), row.names = FALSE)

write.csv(as.data.frame(simple_go_cc_T338I), file.path(dirOUT, "compareCluster_simpleGO_CC_T338I.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_cc_T338M), file.path(dirOUT, "compareCluster_simpleGO_CC_T338M.csv"), row.names = FALSE)
