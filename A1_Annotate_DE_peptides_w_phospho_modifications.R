### Add the new columns in the "combined_modified_peptide.csv" to show differential expression

library(tidyverse)

# The directories

dir250527 = "/media/akira/argentee/proteome/250527_FragPipeAnalystR_1.0.5_MaxLFQ_Proteome_Prashant_Reorganized/20250527"

## Identify the data file that contains PTM info for peptides 
mod_peptide_file = "combined_modified_peptide.tsv"

dir00 <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dirDE = file.path(dir00, "DE_results")
dirK = file.path(dir00, "Kclust_all")

dirCONS = file.path(dir00, "consensus")
if (!dir.exists(dirCONS)) dir.create(dirCONS)
# dirFASTA = file.path(dirCONS, "FASTA_DE")
# if (!dir.exists(dirFASTA)) dir.create(dirFASTA)

dirPAT = file.path(dirCONS, "phospho_patterns_DE")
if (!dir.exists(dirPAT)) dir.create(dirPAT)

# Load peptide data
peptide_df <- read_tsv(file.path(dir250527, mod_peptide_file))
colnames(peptide_df) <- gsub(" ", "_", colnames(peptide_df))

#### Update the obsolete gene symbols TTC26 and C5orf51
peptide_df$Gene[peptide_df$Gene == "TTC26"] = "IFT56"
peptide_df$Gene[peptide_df$Gene == "C5orf51"] = "RIMOC1"

#### Making a small subset df only contains key ID columns without quantitative value columns
peptide_subset <- peptide_df %>%
  select(Peptide_Sequence, Modified_Sequence, Gene, Protein_ID, Protein_Description)

################ Step 1: Annotate the combined_modified_peptide.csv with the DE_bpca_vsn_50_peptides.csv -- the new differential expression output for peptides
### Select key columns and logical columns for DE
DE_peptide_csv = "DE_bpca_vsn_50_result_pep.csv"
de_pep = read.csv(file.path(dirDE, DE_peptide_csv))
de_pep_select = de_pep %>%
  select("Peptide_Sequence",
         "Protein_ID",
         "Gene",
         "Description",
         ends_with("_significant")
         )

# First, identify the "_significant" columns
signif_cols <- grep("_significant$", names(de_pep_select), value = TRUE)

# Then, perform a left join to bring in the "_significant" columns
peptide_tagged <- peptide_subset %>%
  left_join(de_pep_select %>%
              select(Peptide_Sequence, all_of(signif_cols)),
            by = "Peptide_Sequence")


#### Select pY containing sequences in "peptide_select"
library(stringr)

### Peptides with phosphorylated residues
pY_peptides <- peptide_tagged %>%
  filter(str_detect(Modified_Sequence, "Y\\[79\\.9663\\]"))

pT_peptides <- peptide_tagged %>%
  filter(str_detect(Modified_Sequence, "T\\[79\\.9663\\]"))

pS_peptides <- peptide_tagged %>%
  filter(str_detect(Modified_Sequence, "S\\[79\\.9663\\]"))


## Step 1: Summarize each phospho dataset
# IMPORTANT: Handle multiple entries by summarizing with any() to ensure that if any instance of a peptide is phosphorylated on a residue, it is marked TRUE
pY_summary <- pY_peptides %>%
  group_by(Peptide_Sequence) %>%
  summarise(pTyr = any(TRUE), .groups = 'drop')

pT_summary <- pT_peptides %>%
  group_by(Peptide_Sequence) %>%
  summarise(pThr = any(TRUE), .groups = 'drop')

pS_summary <- pS_peptides %>%
  group_by(Peptide_Sequence) %>%
  summarise(pSer = any(TRUE), .groups = 'drop')

## Step 2: Start with your main dataframe
annotated_df <- de_pep_select %>%
  left_join(pY_summary, by = "Peptide_Sequence") %>%
  left_join(pT_summary, by = "Peptide_Sequence") %>%
  left_join(pS_summary, by = "Peptide_Sequence")

## Step 3: Replace NAs with FALSE
annotated_df <- annotated_df %>%
  mutate(across(c(pTyr, pThr, pSer), ~replace_na(.x, FALSE)))


######### Draw a bar graph for phosphorylation co-occurrence patterns #################################

#### Total Phospho Patterns 
# Step 1: Create a combination column
# The Combination column encodes the presence (1) or absence (0) of each phosphorylation type in the order: pTyr, pThr, pSer
annotated_df <- annotated_df %>%
  mutate(Combination = paste0(as.integer(pTyr), as.integer(pThr), as.integer(pSer)))

write.csv(annotated_df, file.path(dirCONS, "DE_bpca_vsn_50_peptide_phosphoMAP.csv"), row.names = F)

# Step 2: Count each combination
combination_counts <- annotated_df %>%
  count(Combination, name = "Count")

# Step 3: Plot the bar chart
pdf(file.path(dirPAT, "Phospho_patterns_total.pdf"))
print(
  ggplot(combination_counts, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "Co-occurrence Patterns of Phosphorylation Types",
       x = "Phosphorylation Combination (pTyr, pThr, pSer)",
       y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
)
dev.off()

#### Comparison Set A: GFP vs P2E8 or GFP vs WT phospho patterns
# Step 1 is the same as Total Phospho Patterns above
# Step 2: Count each combination for DE
combination_counts_P8E2 <- annotated_df %>%
  filter(GFP_vs_P8E2_significant == TRUE) %>%
  count(Combination, name = "Count")

combination_counts_WT <- annotated_df %>%
  filter(GFP_vs_WT_significant == TRUE) %>%
  count(Combination, name = "Count")

# Step 3: Plot the bar chart for DE
pdf(file.path(dirPAT,"Phospho_patterns_P8E2.pdf"))
print(
  ggplot(combination_counts_P8E2, aes(x = Combination, y = Count, fill = Combination)) +
    geom_bar(stat = "identity") +
    labs(title = "Co-occurrence Patterns of Phosphorylation Types",
         x = "Phosphorylation Combination (pTyr, pThr, pSer)",
         y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
)
dev.off()

pdf(file.path(dirPAT,"Phospho_patterns_WT.pdf"))
print(
  ggplot(combination_counts_WT, aes(x = Combination, y = Count, fill = Combination)) +
    geom_bar(stat = "identity") +
    labs(title = "Co-occurrence Patterns of Phosphorylation Types",
         x = "Phosphorylation Combination (pTyr, pThr, pSer)",
         y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
)
dev.off()

# Step 4: a modified Step 3 to put two plots side-by-side
library(patchwork)

# Define a consistent color palette for all combinations
all_combinations <- union(combination_counts_WT$Combination, combination_counts_P8E2$Combination)
color_palette <- RColorBrewer::brewer.pal(n = length(all_combinations), name = "Set2")
names(color_palette) <- all_combinations

# Create the plots with manual fill scale
plot_P8E2 <- ggplot(combination_counts_P8E2, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "GFP vs P8E2", x = "Phosphorylation Pattern", y = "Peptides") +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  ylim(0, max(c(combination_counts_P8E2$Count, combination_counts_WT$Count)))

plot_WT <- ggplot(combination_counts_WT, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "GFP vs WT", x = "Phosphorylation Pattern", y = "Peptides") +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  ylim(0, max(c(combination_counts_P8E2$Count, combination_counts_WT$Count)))

# Combine and save
pdf(file.path(dirPAT, "Phospho_patterns_combined_P8E2_WT.pdf"), width = 12, height = 6)
print(plot_P8E2 + plot_WT + plot_layout(ncol = 2))
dev.off()

#### Comparison Set B: GFP vs T338I or GFP vs T338M phospho patterns
# Step 1 is the same as Total Phospho Patterns above
# Step 2: Count each combination for DE
combination_counts_T338I <- annotated_df %>%
  filter(GFP_vs_WT_T338I_significant == TRUE) %>%
  count(Combination, name = "Count")

combination_counts_T338M <- annotated_df %>%
  filter(GFP_vs_WT_T338M_significant == TRUE) %>%
  count(Combination, name = "Count")

# Step 3: Plot the bar chart for DE
pdf(file.path(dirPAT,"Phospho_patterns_T338I.pdf"))
print(
  ggplot(combination_counts_T338I, aes(x = Combination, y = Count, fill = Combination)) +
    geom_bar(stat = "identity") +
    labs(title = "Co-occurrence Patterns of Phosphorylation Types",
         x = "Phosphorylation Combination (pTyr, pThr, pSer)",
         y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
)
dev.off()

pdf(file.path(dirPAT,"Phospho_patterns_T338M.pdf"))
print(
  ggplot(combination_counts_T338M, aes(x = Combination, y = Count, fill = Combination)) +
    geom_bar(stat = "identity") +
    labs(title = "Co-occurrence Patterns of Phosphorylation Types",
         x = "Phosphorylation Combination (pTyr, pThr, pSer)",
         y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
)
dev.off()

# Step 4: a modified Step 3 to put two plots side-by-side
# Define a consistent color palette for all combinations
all_combinations_T <- union(combination_counts_T338I$Combination, combination_counts_T338M$Combination)
color_palette <- RColorBrewer::brewer.pal(n = length(all_combinations_T), name = "Set2")
names(color_palette) <- all_combinations_T

# Create the plots with manual fill scale
plot_T338I <- ggplot(combination_counts_T338I, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "WT_T338I", x = "Phosphorylation Pattern", y = "Peptides") +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  ylim(0, max(c(combination_counts_T338I$Count, combination_counts_T338M$Count)))

plot_T338M <- ggplot(combination_counts_T338M, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "WT_T338M", x = "Phosphorylation Pattern", y = "Peptides") +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  ylim(0, max(c(combination_counts_T338I$Count, combination_counts_T338M$Count)))

# Combine and save
pdf(file.path(dirPAT, "Phospho_patterns_combined_T338I_T338M.pdf"), width = 12, height = 6)
print(plot_T338I + plot_T338M + plot_layout(ncol = 2))
dev.off()

#### Comparison Set C: GFP vs P8E2_T338I or GFP vs P8E2_T338M phospho patterns
# Step 1 is the same as Total Phospho Patterns above
# Step 2: Count each combination for DE
combination_counts_P8E2xT338I <- annotated_df %>%
  filter(GFP_vs_P8E2_T338I_significant == TRUE) %>%
  count(Combination, name = "Count")

combination_counts_P8E2xT338M <- annotated_df %>%
  filter(GFP_vs_P8E2_T338M_significant == TRUE) %>%
  count(Combination, name = "Count")

# Step 3: Plot the bar chart for DE
pdf(file.path(dirPAT,"Phospho_patterns_P8E2_T338I.pdf"))
print(
  ggplot(combination_counts_P8E2xT338I, aes(x = Combination, y = Count, fill = Combination)) +
    geom_bar(stat = "identity") +
    labs(title = "Co-occurrence Patterns of Phosphorylation Types",
         x = "Phosphorylation Combination (pTyr, pThr, pSer)",
         y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
)
dev.off()

pdf(file.path(dirPAT,"Phospho_patterns_P8E2_T338M.pdf"))
print(
  ggplot(combination_counts_P8E2xT338M, aes(x = Combination, y = Count, fill = Combination)) +
    geom_bar(stat = "identity") +
    labs(title = "Co-occurrence Patterns of Phosphorylation Types",
         x = "Phosphorylation Combination (pTyr, pThr, pSer)",
         y = "Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2")
)
dev.off()

# Step 4: a modified Step 3 to put two plots side-by-side
# Define a consistent color palette for all combinations
all_combinations_P <- union(combination_counts_P8E2xT338I$Combination, combination_counts_P8E2xT338M$Combination)
color_palette <- RColorBrewer::brewer.pal(n = length(all_combinations_P), name = "Set2")
names(color_palette) <- all_combinations_P

# Create the plots with manual fill scale
plot_P8E2xT338I <- ggplot(combination_counts_P8E2xT338I, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "P8E2_T338I", x = "Phosphorylation Pattern", y = "Peptides") +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  ylim(0, max(c(combination_counts_P8E2xT338I$Count, combination_counts_P8E2xT338M$Count)))

plot_P8E2xT338M <- ggplot(combination_counts_P8E2xT338M, aes(x = Combination, y = Count, fill = Combination)) +
  geom_bar(stat = "identity") +
  labs(title = "P8E2_T338M", x = "Phosphorylation Pattern", y = "Peptides") +
  theme_minimal() +
  scale_fill_manual(values = color_palette) +
  ylim(0, max(c(combination_counts_P8E2xT338I$Count, combination_counts_P8E2xT338M$Count)))

# Combine and save
pdf(file.path(dirPAT, "Phospho_patterns_combined_P8E2xT338I_P8E2xT338M.pdf"), width = 12, height = 6)
print(plot_P8E2xT338I + plot_P8E2xT338M + plot_layout(ncol = 2))
dev.off()

