#######################################################################
# SET OPERATION for Venn's Diagram for Differential Expression Analysis
#######################################################################

# PARTs 1-2 --> RUN to generate csv files as input files for the custom Python script wrappers v2p and v3p to generate Venn's diagrams. RUN PART 1 and PART 2. To export the count data for additional comparison sets, add them as needed to PART 2, Step 5

# PART 3 (starting ~ line 270) --> Optional Post-Hoc subset analysis to identify the proteins/peptides in subsets. RUN PART 1, PART 2 (STEP 1 and STEP2 only) then PART 3

###### PART 1: Before SET OPERATIONS ##########################################
# IMPORTANT --> CUSTOMIZE for your data structure and DEFINE your comparison sets

### Identify the files for data
dir0 = "/media/akira/argentee/proteome/251014_FragPipeAnalystR_1.1.0_MaxLFQ_250527"
dir1 = file.path(dir0, "DE_results")
pro = "DE_bpca_vsn_50_result_pro.csv" 
pep = "DE_bpca_vsn_50_result_pep.csv"

df_pro = read.csv(file.path(dir1, pro))
df_pep = read.csv(file.path(dir1, pep))

### Define the sets for Venn's diagram
set_list1 = c("GFP_vs_WT", "GFP_vs_P8E2")
set_list2 = c("GFP_vs_WT", "GFP_vs_WT_T338I", "GFP_vs_WT_T338M")
set_list3 = c("GFP_vs_P8E2", "GFP_vs_P8E2_T338I", "GFP_vs_P8E2_T338M")

### PART 2: Find DE proteins or peptides in the set lists #######################
## First select proteins or peptides that are "TRUE" in the corresponding columns named "*_significant"

#### STEP 1: Define Functions for Set Creation

## Function to get ALL DE items (combined) based only on significance
get_all_de_sets <- function(data, comparison_list, id_column) {
  de_sets <- list()
  for (comparison in comparison_list) {
    sig_col_name <- paste0(comparison, "_significant")
    if (sig_col_name %in% names(data)) {
      is_significant <- !is.na(data[[sig_col_name]]) & data[[sig_col_name]]
      ids <- data[is_significant, id_column]
      de_sets[[comparison]] <- ids
    } else {
      warning(paste("Column '", sig_col_name, "' not found.", sep=""))
    }
  }
  return(de_sets)
}

## Function to get DIRECTIONAL DE items (Upregulated or Downregulated)
get_directional_sets <- function(data, comparison_list, id_column, direction) {
  de_sets <- list()
  for (comparison in comparison_list) {
    sig_col_name <- paste0(comparison, "_significant")
    diff_col_name <- paste0(comparison, "_diff")
    if (sig_col_name %in% names(data) && diff_col_name %in% names(data)) {
      is_significant <- !is.na(data[[sig_col_name]]) & data[[sig_col_name]]
      significant_subset <- data[is_significant, ]
      
      ids <- c()
      if (direction == "Upregulated") {
        is_selected <- !is.na(significant_subset[[diff_col_name]]) & significant_subset[[diff_col_name]] < 0
        ids <- significant_subset[is_selected, id_column]
      } else if (direction == "Downregulated") {
        is_selected <- !is.na(significant_subset[[diff_col_name]]) & significant_subset[[diff_col_name]] > 0
        ids <- significant_subset[is_selected, id_column]
      }
      de_sets[[comparison]] <- ids
    } else {
      if (!sig_col_name %in% names(data)) warning(paste("Column '", sig_col_name, "' not found.", sep=""))
      if (!diff_col_name %in% names(data)) warning(paste("Column '", diff_col_name, "' not found.", sep=""))
    }
  }
  return(de_sets)
}

#### STEP 2: Apply Functions to Create Sets

## First, create the unique peptide identifier
df_pep$Gene_Peptide_ID <- paste(df_pep$Gene, df_pep$Peptide_Sequence, sep = "_")

## --- NEW SECTION: Get COMBINED Sets (Up & Down Together) ---
## For Proteins (3 combined sets)
protein_sets1_all <- get_all_de_sets(df_pro, set_list1, "Gene")
protein_sets2_all <- get_all_de_sets(df_pro, set_list2, "Gene")
protein_sets3_all <- get_all_de_sets(df_pro, set_list3, "Gene")

## For Peptides (3 combined sets)
peptide_sets1_all <- get_all_de_sets(df_pep, set_list1, "Gene_Peptide_ID")
peptide_sets2_all <- get_all_de_sets(df_pep, set_list2, "Gene_Peptide_ID")
peptide_sets3_all <- get_all_de_sets(df_pep, set_list3, "Gene_Peptide_ID")


## --- Get SEPARATE UP/DOWN Sets ---
## For Proteins (directional sets)
protein_sets1_up <- get_directional_sets(df_pro, set_list1, "Gene", "Upregulated")
protein_sets1_down <- get_directional_sets(df_pro, set_list1, "Gene", "Downregulated")
protein_sets2_up <- get_directional_sets(df_pro, set_list2, "Gene", "Upregulated")
protein_sets2_down <- get_directional_sets(df_pro, set_list2, "Gene", "Downregulated")
protein_sets3_up <- get_directional_sets(df_pro, set_list3, "Gene", "Upregulated")
protein_sets3_down <- get_directional_sets(df_pro, set_list3, "Gene", "Downregulated")

## For Peptides (directional sets)
peptide_sets1_up <- get_directional_sets(df_pep, set_list1, "Gene_Peptide_ID", "Upregulated")
peptide_sets1_down <- get_directional_sets(df_pep, set_list1, "Gene_Peptide_ID", "Downregulated")
peptide_sets2_up <- get_directional_sets(df_pep, set_list2, "Gene_Peptide_ID", "Upregulated")
peptide_sets2_down <- get_directional_sets(df_pep, set_list2, "Gene_Peptide_ID", "Downregulated")
peptide_sets3_up <- get_directional_sets(df_pep, set_list3, "Gene_Peptide_ID", "Upregulated")
peptide_sets3_down <- get_directional_sets(df_pep, set_list3, "Gene_Peptide_ID", "Downregulated")

#### STEP 3: Verification (Optional)
## Check the new 'all' sets
# print("Number of ALL DE proteins in set 1 (Up & Down combined):")
# print(lapply(protein_sets1_all, length))
# 
# ## Check the directional sets (as before)
# print("Number of Upregulated proteins in set 1:")
# print(lapply(protein_sets1_up, length))


#### STEP 4: Perform Set Operations and Export Counts for Venn Diagrams

# First, let's create a dedicated directory for the output CSV files.
dir_out = file.path(dir1, "venn_counts")
dir.create(dir_out, showWarnings = FALSE)

## --------------------------------------------------------------------------
## Function 1: Calculate counts for a 2-set comparison
## --------------------------------------------------------------------------
calculate_venn2_counts <- function(sets_list) {
  # Check if the input is a list with exactly two sets
  if (!is.list(sets_list) || length(sets_list) != 2) {
    stop("Input must be a list containing exactly two sets (vectors).")
  }
  
  # Assign sets to simple names for clarity
  A <- sets_list[[1]]
  B <- sets_list[[2]]
  
  # Perform set operations
  A_only <- setdiff(A, B) # Corresponds to 'Ab'
  B_only <- setdiff(B, A) # Corresponds to 'aB'
  A_and_B <- intersect(A, B) # Corresponds to 'AB'
  
  # Get the counts for each segment
  counts <- c(
    length(A_only),
    length(B_only),
    length(A_and_B)
  )
  
  # --- MODIFIED PART ---
  # Name the counts using the specified 'Ab' convention
  names(counts) <- c("Ab", "aB", "AB")
  
  return(as.data.frame(t(counts)))
}

## --------------------------------------------------------------------------
## Function 2: Calculate counts for a 3-set comparison
## --------------------------------------------------------------------------
calculate_venn3_counts <- function(sets_list) {
  # Check if the input is a list with exactly three sets
  if (!is.list(sets_list) || length(sets_list) != 3) {
    stop("Input must be a list containing exactly three sets (vectors).")
  }
  
  # Assign sets to simple names
  A <- sets_list[[1]]
  B <- sets_list[[2]]
  C <- sets_list[[3]]
  
  # Perform set operations to find the 7 unique segments
  A_only <- setdiff(A, union(B, C)) # Corresponds to 'Abc'
  B_only <- setdiff(B, union(A, C)) # Corresponds to 'aBc'
  C_only <- setdiff(C, union(A, B)) # Corresponds to 'abC'
  A_B_only <- setdiff(intersect(A, B), C) # Corresponds to 'ABc'
  A_C_only <- setdiff(intersect(A, C), B) # Corresponds to 'AbC'
  B_C_only <- setdiff(intersect(B, C), A) # Corresponds to 'aBC'
  A_B_C <- intersect(intersect(A, B), C) # Corresponds to 'ABC'
  
  # Get counts for each segment, ordered as in your example
  counts <- c(
    length(A_only),
    length(B_only),
    length(A_B_only), # Note the order change to match your example (ABc before abC)
    length(C_only),
    length(A_C_only),
    length(B_C_only),
    length(A_B_C)
  )
  
  # --- MODIFIED PART ---
  # Name the counts using the specified 'Abc' convention and order
  names(counts) <- c("Abc", "aBc", "ABc", "abC", "AbC", "aBC", "ABC")
  
  return(as.data.frame(t(counts)))
}

## --------------------------------------------------------------------------
## Main Wrapper Function: Analyze a list and export the CSV
## --------------------------------------------------------------------------
export_venn_counts <- function(sets_list, filename) {
  # Determine if it's a 2-set or 3-set comparison
  num_sets <- length(sets_list)
  
  if (num_sets == 2) {
    cat("Performing 2-set analysis on:", deparse(substitute(sets_list)), "\n")
    counts_df <- calculate_venn2_counts(sets_list)
  } else if (num_sets == 3) {
    cat("Performing 3-set analysis on:", deparse(substitute(sets_list)), "\n")
    counts_df <- calculate_venn3_counts(sets_list)
  } else {
    warning(paste("Skipping:", deparse(substitute(sets_list)), "- it does not contain 2 or 3 sets."))
    return(NULL)
  }
  
  # Construct the full file path and write the CSV
  full_path <- file.path(dir_out, filename)
  write.csv(counts_df, full_path, row.names = FALSE)
  cat("Results saved to:", full_path, "\n\n")
}


#### STEP 5: Apply the Functions to Your Datasets

## Without separating the directions
# Proteins
export_venn_counts(
  sets_list = protein_sets1_all, 
  filename = "Venn2_protein_GFPvWT_GFPvP8E2_all.csv"
)

export_venn_counts(
  sets_list = protein_sets2_all, 
  filename = "Venn3_protein_GFPvWT_GFPvT338I_GFPvT338M_all.csv"
)

export_venn_counts(
  sets_list = protein_sets3_all, 
  filename = "Venn3_protein_GFPvP8_GFPvP8T338I_GFPvP8T338M_all.csv"
)

# Peptides
export_venn_counts(
  sets_list = peptide_sets1_all, 
  filename = "Venn2_peptide_GFPvWT_GFPvP8E2_all.csv"
)

export_venn_counts(
  sets_list = peptide_sets2_all, 
  filename = "Venn3_peptide_GFPvWT_GFPvT338I_GFPvT338M_all.csv"
)

export_venn_counts(
  sets_list = peptide_sets3_all, 
  filename = "Venn3_peptide_GFPvP8_GFPvP8T338I_GFPvP8T338M_all.csv"
)

## Separating upregulated or downregulated
# protein_sets1_up 
export_venn_counts(
  sets_list = protein_sets1_up, 
  filename = "Venn2_protein_GFPvWT_GFPvP8E2_upregulated.csv"
)

# protein_sets1_down
export_venn_counts(
  sets_list = protein_sets1_down, 
  filename = "Venn2_protein_GFPvWT_GFPvP8E2_downregulated.csv"
)

##############################################################################
# PART 3: (OPTIONAL) POST-HOC SUBSET ANALYSIS
##############################################################################

# Create a new directory for the subset list files
dir_subsets = file.path(dir1, "venn_subsets")
dir.create(dir_subsets, showWarnings = FALSE)

## --------------------------------------------------------------------------
## Function 1: Get subset lists for a 2-set comparison
## --------------------------------------------------------------------------
get_venn2_subsets <- function(sets_list) {
  if (!is.list(sets_list) || length(sets_list) != 2) {
    stop("Input must be a list containing exactly two sets.")
  }
  
  # Assign sets to simple names
  A <- sets_list[[1]]
  B <- sets_list[[2]]
  
  # Perform set operations to get the identifiers in each subset
  Ab <- setdiff(A, B)
  aB <- setdiff(B, A)
  AB <- intersect(A, B)
  
  # Return a named list of the actual identifier vectors
  subsets <- list(
    Ab = Ab,
    aB = aB,
    AB = AB
  )
  return(subsets)
}

## --------------------------------------------------------------------------
## Function 2: Get subset lists for a 3-set comparison
## --------------------------------------------------------------------------
get_venn3_subsets <- function(sets_list) {
  if (!is.list(sets_list) || length(sets_list) != 3) {
    stop("Input must be a list containing exactly three sets.")
  }
  
  # Assign sets to simple names
  A <- sets_list[[1]]
  B <- sets_list[[2]]
  C <- sets_list[[3]]
  
  # Perform set operations to find the 7 unique segments
  Abc <- setdiff(A, union(B, C))
  aBc <- setdiff(B, union(A, C))
  abC <- setdiff(C, union(A, B))
  ABc <- setdiff(intersect(A, B), C)
  AbC <- setdiff(intersect(A, C), B)
  aBC <- setdiff(intersect(B, C), A)
  ABC <- intersect(intersect(A, B), C)
  
  # Return a named list of the identifier vectors in the desired order
  subsets <- list(
    Abc = Abc,
    aBc = aBc,
    ABc = ABc,
    abC = abC,
    AbC = AbC,
    aBC = aBC,
    ABC = ABC
  )
  return(subsets)
}

## --------------------------------------------------------------------------
## Main Wrapper Function: Get subsets and write to a formatted text file
## --------------------------------------------------------------------------
export_venn_subsets <- function(sets_list, filename) {
  num_sets <- length(sets_list)
  
  if (num_sets == 2) {
    cat("Generating 2-set subsets for:", deparse(substitute(sets_list)), "\n")
    subsets <- get_venn2_subsets(sets_list)
  } else if (num_sets == 3) {
    cat("Generating 3-set subsets for:", deparse(substitute(sets_list)), "\n")
    subsets <- get_venn3_subsets(sets_list)
  } else {
    warning(paste("Skipping:", deparse(substitute(sets_list)), "- does not contain 2 or 3 sets."))
    return(NULL)
  }
  
  # Write the results to a single, organized text file
  full_path <- file.path(dir_subsets, filename)
  sink(full_path) # Redirects console output to the file
  
  for (subset_name in names(subsets)) {
    cat(paste("## Subset:", subset_name, "(Count:", length(subsets[[subset_name]]), ")\n"))
    # Write the identifiers, one per line
    writeLines(subsets[[subset_name]])
    cat("\n") # Add a blank line for readability
  }
  
  sink() # Closes the connection, redirecting output back to the console
  cat("Subset lists saved to:", full_path, "\n\n")
}

#### STEP 7: Apply the Functions to Your Datasets

# Prerequisite: Ensure the set lists from the previous script are available in your R environment
# (e.g., protein_sets1_all, protein_sets2_up, etc.)

## --- protein_sets1_all (GFPvWT, GFPvP8) ---
export_venn_subsets(
  sets_list = protein_sets1_all, 
  filename = "Subsets_protein_all_GFPvWT_GFPvP8.txt"
)

## --- peptide_sets1_all (GFPvWT, GFPvP8) ---
export_venn_subsets(
  sets_list = peptide_sets1_all, 
  filename = "Subsets_peptide_all_GFPvWT_GFPvP8.txt"
)

## protein_sets2_all ()
export_venn_subsets(
  sets_list = protein_sets2_all, 
  filename = "Subsets_protein_all_GFPvWT_GFPvT338I_GFPvT338M.txt"
)

## protein_sets3_all ()
export_venn_subsets(
  sets_list = protein_sets3_all, 
  filename = "Subsets_protein_all_GFPvP8_GFPvP8T338I_GFPvP8T338M.txt"
)

##### OPTIONAL
## --- Example for a 3-set list (protein_sets2_up) ---
export_venn_subsets(
  sets_list = protein_sets2_up, 
  filename = "Subsets_protein_upregulated_set2.txt"
)

## --- Example for another 3-set list (peptide_sets3_down) ---
export_venn_subsets(
  sets_list = peptide_sets3_down,
  filename = "Subsets_peptide_downregulated_set3.txt"
)
