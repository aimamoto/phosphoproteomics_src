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

############### The Directories and Data
dir0 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dirK <- file.path(dir0, "Kclust_all")

clustDA <- "clusters_bpca_vsn_50_result_pro.csv"

dfclust <- read.csv(file.path(dirK, clustDA))

dirOUT <- file.path(dirK, "enrich_by_kClusters")
if (!dir.exists(dirOUT)) dir.create(dirOUT)

############### Map UniProt IDs to Entrez and Symbol
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


# ---- Prepare Gene List for compareCluster ----
gene_list <- dfclust_mapped %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(Cluster) %>%
  summarise(genes = list(unique(ENTREZID))) %>%
  deframe()

# ---- Run compareCluster for Multiple Ontologies ----

# GO: combined BP, MF, and CC
# compare_go_all = compareCluster(geneCluster = gene_list,
#                                 fun = "enrichGO",
#                                 OrgDb = org.Hs.eg.db,
#                                 ont = "ALL",
#                                 pvalueCutoff = 0.01)
# 
#     #### Modify compare_go_all to apply a filter by GeneRatio > 0.05
#       # Step 1: Extract and filter
#         filtered_df <- compare_go_all@compareClusterResult %>%
#           separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE) %>%
#           mutate(GeneRatio_numeric = numerator / denominator) %>%
#           filter(GeneRatio_numeric > 0.05) %>%
#           # Restore GeneRatio
#           mutate(GeneRatio = paste(numerator, denominator, sep = "/")) %>% 
#           select(-numerator, -denominator, -GeneRatio_numeric)
# 
#       # Step 2: Replace the slot
#       compare_go_all_filtered <- compare_go_all
#       compare_go_all_filtered@compareClusterResult <- filtered_df

# GO-BP
compare_go_bp <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "BP",
                                pvalueCutoff = 0.01)
simple_go_bp = simplify(compare_go_bp)

# GO-MF
compare_go_mf <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "MF",
                                pvalueCutoff = 0.01)
simple_go_mf = simplify(compare_go_mf)

# GO-CC
compare_go_cc <- compareCluster(geneCluster = gene_list,
                                fun = "enrichGO",
                                OrgDb = org.Hs.eg.db,
                                ont = "CC",
                                pvalueCutoff = 0.01)
simple_go_cc = simplify(compare_go_cc)

# KEGG
compare_kegg <- compareCluster(geneCluster = gene_list,
                               fun = "enrichKEGG",
                               organism = "hsa",
                               pvalueCutoff = 0.01)

# MSigDB Hallmark
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
compare_hallmark <- compareCluster(geneCluster = gene_list,
                                   fun = "enricher",
                                   TERM2GENE = msig_hallmark %>% select(gs_name, ncbi_gene),
                                   pvalueCutoff = 0.01)

# WikiPathways
msig_wiki <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
compare_wiki <- compareCluster(geneCluster = gene_list,
                               fun = "enricher",
                               TERM2GENE = msig_wiki %>% select(gs_name, ncbi_gene),
                               pvalueCutoff = 0.01)

############ ---- Save Dotplots ----

pdf(file.path(dirOUT, "compareCluster_KEGG.pdf"), width = 10, height = 10)
print(
  dotplot(compare_kegg, showCategory = 10) + 
    ggtitle("KEGG by Cluster for All Groups Combined") +
    theme(
      axis.text.y = element_text(size = 15), 
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

pdf(file.path(dirOUT, "compareCluster_Hallmark.pdf"), width = 10, height = 8)
print(
  dotplot(compare_hallmark, showCategory = 10) + 
    ggtitle("MsigDB Hallmark by Cluster for All Groups Combined") +
    theme(
      axis.text.y = element_text(size = 15), 
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

pdf(file.path(dirOUT, "compareCluster_WikiPathways.pdf"), width = 10, height = 8)
print(
  dotplot(compare_wiki, showCategory = 10) + 
    ggtitle("WikiPathways by Cluster for All Groups Combined") +
    theme(
      axis.text.y = element_text(size = 15), 
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

# simplified GO
pdf(file.path(dirOUT, "compareCluster_simpleGO_BP.pdf"), width = 10, height = 15)
print(dotplot(simple_go_bp, showCategory = 8) + 
        ggtitle("Simplified GO-BP by Cluster for All Groups Combined") + 
        theme(
          axis.text.y = element_text(size = 15), 
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

pdf(file.path(dirOUT, "compareCluster_simpleGO_CC.pdf"), width = 10, height = 15)
print(dotplot(simple_go_cc, showCategory = 10) + 
        ggtitle("Simplified GO-CC by Cluster for All Groups Combined") + 
        theme(
          axis.text.y = element_text(size = 15), 
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

pdf(file.path(dirOUT, "compareCluster_simpleGO_MF.pdf"), width = 10, height = 15)
print(dotplot(simple_go_mf, showCategory = 10) + 
        ggtitle("Simplified GO Molecular Function") +
        theme(
          axis.text.y = element_text(size = 15), 
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

# ---- Save Tables ----
# write.csv(as.data.frame(compare_go_all), file.path(dirOUT, "compareCluster_GO_ALL.csv"), row.names = FALSE)

# write.csv(as.data.frame(compare_go_bp), file.path(dirOUT, "compareCluster_GO_BP.csv"), row.names = FALSE)
# write.csv(as.data.frame(compare_go_mf), file.path(dirOUT, "compareCluster_GO_MF.csv"), row.names = FALSE)
# write.csv(as.data.frame(compare_go_cc), file.path(dirOUT, "compareCluster_GO_CC.csv"), row.names = FALSE)
# write.csv(as.data.frame(compare_kegg), file.path(dirOUT, "compareCluster_KEGG.csv"), row.names = FALSE)

write.csv(as.data.frame(compare_hallmark), file.path(dirOUT, "compareCluster_Hallmark.csv"), row.names = FALSE)
write.csv(as.data.frame(compare_wiki), file.path(dirOUT, "compareCluster_WikiPathways.csv"), row.names = FALSE)

write.csv(as.data.frame(simple_go_bp), file.path(dirOUT, "compareCluster_simpleGO_BP.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_mf), file.path(dirOUT, "compareCluster_simpleGO_MF.csv"), row.names = FALSE)
write.csv(as.data.frame(simple_go_cc), file.path(dirOUT, "compareCluster_simpleGO_CC.csv"), row.names = FALSE)