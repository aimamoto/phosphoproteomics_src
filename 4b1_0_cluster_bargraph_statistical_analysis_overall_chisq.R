# --- Load and filter data as before ---
dirDATA <- "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527/Kclust_all/kCluster_bargraph"
data <- "Final_kCluster_Counts.csv"
df <- read.csv(file.path(dirDATA, data))
df_filtered <- df[df$Cluster != "NoData", ]

# --- Create the main contingency tables ---
peptide_table <- xtabs(Count ~ Group + Cluster, 
                       data = subset(df_filtered, Type == "peptides"))

protein_table <- xtabs(Count ~ Group + Cluster, 
                       data = subset(df_filtered, Type == "proteins"))

# --- Run the Overall Fisher's Exact Test ---
# Fisher's test is recommended over Chi-squared due to potentially sparse counts.
# Note: This may take a moment to compute for larger tables.
# --- Run the Overall Test with Simulation ---

cat("--- Overall Test for Peptide Data (with simulation) ---\n")
# The only change is adding simulate.p.value = TRUE
overall_peptide_test <- chisq.test(peptide_table, simulate.p.value = TRUE, B = 10000)
print(overall_peptide_test)

cat("\n--- Overall Test for Protein Data (with simulation) ---\n")
# You may need to add it here as well if the protein table is also large
overall_protein_test <- chisq.test(protein_table, simulate.p.value = TRUE, B = 10000)
print(overall_protein_test)