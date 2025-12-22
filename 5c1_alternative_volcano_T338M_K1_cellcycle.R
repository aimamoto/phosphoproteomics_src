# Make sure libraries are loaded
library(ggrepel)
library(FragPipeAnalystR)
library(SummarizedExperiment)
library(tidyverse)

#### Load data
dir0 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
se_save_dir <- file.path(dir0, "DE_results")
load(file.path(se_save_dir, "de_results_bpca_list.RData"))

#### Load and filter the gene list for annotation
dirK = file.path(dir0, "Kclust_all")
dirKclustGO <- file.path(dirK, "GO_geneSymbols_enrich_by_kClusters")
BP_T338M_symbol <- "compareCluster_simpleGO_BP_T338M_by_Symbol.csv"

BP_T338M_df <- read.csv(file.path(dirKclustGO, BP_T338M_symbol))

K1_cellcycle_T338M <- BP_T338M_df %>%
  filter(Cluster == 1) %>%
  filter(str_detect(Description, "mitotic cell cycle"))

################# Focused Volcano Plot Generation
# Define the specific result name and contrast we want to plot
result_name <- "bpca_vsn_50_result_pro"
contrast_name <- "GFP_vs_WT_T338M"

dirOUT = file.path(dirK, "volcano_focused_cellcycle_K1")
if(!dir.exists(dirOUT)) dir.create(dirOUT)

# 1. Directly select the protein-level SE object from the list
de_result <- de_results_bpca_list[[result_name]]

# 2. Extract the full results data from the SE object
results_df <- as.data.frame(rowData(de_result))

# 3. Get the list of gene symbols to label
proteins_to_label <- K1_cellcycle_T338M$geneSymbol

# 4. Create the specific subset of data for labeling
label_data <- results_df %>%
  filter(Gene %in% proteins_to_label)

# Define the output file path
filename <- paste0("volcano_", result_name, "_", contrast_name, "_CellCycleLabeled.pdf")
filepath <- file.path(dirOUT, filename)

# 5. Generate and save the plot
pdf(filepath)

# Create the base plot object
base_plot <- plot_volcano(de_result, contrast_name, name_col = "Gene") +
  ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey40")

# Add the annotation layer with custom color with the x and y parameter scales specified for GFP vs WT_T338M (fixed by copilot)
annotated_plot <- base_plot +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(
      x = GFP_vs_WT_T338M_diff,
      y = -log10(GFP_vs_WT_T338M_p.adj),
      label = Gene
    ),
    color = "navyblue",
    fontface = "bold",
    size = 3, ## Adjust the font size as needed!
    box.padding = 0.5,
    max.overlaps = Inf,
    segment.color = 'grey50'
  )

print(annotated_plot)
dev.off()

# Confirmation message
message(paste("Volcano plot saved to:", filepath))