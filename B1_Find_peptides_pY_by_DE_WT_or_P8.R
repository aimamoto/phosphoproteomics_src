

library(tidyverse)

dir0 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dir1 = file.path(dir0, "consensus")
dirFA = file.path(dir1, "FASTA_DE")
dir.create(dirFA, showWarnings = F)

phosphoMAPfile = "DE_bpca_vsn_50_peptide_phosphoMAP.csv"
df_DE = read.csv(file.path(dir1, phosphoMAPfile))

#### This is a modified B1 script in which I intend to find all peptides with pY in SET WT, SET P8, and SET WTnotP8

df_ALLpY = df_DE %>%
  filter(pTyr == "TRUE")

df_ALLpYWT = df_ALLpY %>%
  filter(GFP_vs_WT_significant == "TRUE")

df_ALLpYP8 = df_ALLpY %>%
  filter(GFP_vs_P8E2_significant == "TRUE")

df_ALLpYWTnotP8 = dplyr::setdiff(df_ALLpYWT, df_ALLpYP8) # alternatively, 

#### Make fasta file for "Peptide_Sequence" in each with header of > "Gene" | "Protein_ID" | "Description. 
# Function to write a dataframe to a FASTA file
write_fasta <- function(df, file_name) {
  lines <- apply(df, 1, function(row) {
    header <- paste0(">", row["Protein_ID"], " | ", row["Gene"], " | ", row["Description"])
    sequence <- row["Peptide_Sequence"]
    c(header, sequence)
  })
  writeLines(unlist(lines), con = file_name)
}

# Selected dataframes for running the function above
pY_sets_selected <- list(
  df_ALLpYWT = df_ALLpYWT,
  df_ALLpYP8 = df_ALLpYP8,
  df_ALLpYWTnotP8 = df_ALLpYWTnotP8
)

# write one fasta per data.frame; file names will be e.g. "ALLpYWT.fasta"
for (nm in names(pY_sets_selected)) {
  df <- pY_sets_selected[[nm]]
  # remove leading "df_" to get the distinct part; change this if you prefer a different naming scheme
  base <- gsub("^df_", "", nm)
  out_file <- file.path(dirFA, paste0(base, ".fasta"))
  write_fasta(df, out_file)
  message("Wrote ", out_file, " (", nrow(df[!is.na(df$Peptide_Sequence) & df$Peptide_Sequence != "", ]), " sequences)")
}
