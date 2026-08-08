################################################################################
# Project : "MicrobExplorer"
# Script  : "Analysis : Phyloseq on standardized counts"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

suppressPackageStartupMessages({
  # Libraries CRAN
  library(data.table)
  library(readxl)
  # Libraries Bioconductor
  library(phyloseq)
})
# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================
source("workflow/scripts/utils/utils_io.R")

# Inputs
DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])[1]
TAXONOMY <- as.character(snakemake@input[["taxonomy"]])[1]

# Outputs
RDS      <- as.character(snakemake@output[["rds"]])[1]

# Parameters
STAND_COL <- tolower(as.character(snakemake@params[["stand_col"]]))[1]

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ==========================================================================
# 1. Automatic Data Type Detection (Reads vs Contigs vs KEGG)
# ==========================================================================
if (grepl("kegg", SOURCE, ignore.case = TRUE)) {
  target_id <- "kegg_id"
  tax_cols  <- c("ec_number", "level_1", "level_2", "level_3", "gene_description")
} else if (grepl("reads", SOURCE, ignore.case = TRUE)) {
  target_id <- "read_id" 
  tax_cols  <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
} else if (grepl("contigs", SOURCE, ignore.case = TRUE)) {
  target_id <- "contig_id"
  tax_cols  <- c("domain", "phylum", "class", "order", "family", "genus", "species")
} else {
  stop(sprintf("❌ Error: Unknown SOURCE value [%s]. Expected 'kegg', 'reads', or 'contigs'.", SOURCE))
}

# ==========================================================================
# 2. Data loading
# ==========================================================================
if (is.null(STAND_COL) || is.na(STAND_COL) || STAND_COL == "") {
  stop("❌ ERROR: 'stand_col' parameter is missing or empty in Snakemake config.")
}

meta_dt <- load_metadata(METADATA)

all_data <- load_tsv_dir_dynamic(
  paths       = DATA, 
  meta_dt     = meta_dt, 
  select_cols = c(target_id, STAND_COL)
)

# Format IDs as character and clean missing values
all_data[, (target_id) := as.character(get(target_id))]
all_data <- all_data[!is.na(get(target_id)) & get(target_id) != "" & get(target_id) != "NA"]

# Ensure target abundance column is numeric prior to aggregation
all_data[, (STAND_COL) := as.numeric(get(STAND_COL))]
# ==========================================================================
# 3. Build Phyloseq Components
# ==========================================================================

# --- A. OTU TABLE ---
formula_str <- as.formula(paste(target_id, "~ sample_id"))

# Safeguard: fun.aggregate = sum handles multiple identical KOs per sample perfectly
otu_dt <- dcast(
  all_data,
  formula_str,
  value.var = STAND_COL,
  fun.aggregate = sum,
  fill = 0
)

otu_mat <- as.matrix(otu_dt, rownames = target_id)
mode(otu_mat) <- "numeric"

n_na <- sum(is.na(otu_mat))
if (n_na > 0) {
  warning(sprintf("⚠️ WARNING: %d NA values detected in OTU matrix — replaced with 0.", n_na))
  otu_mat[is.na(otu_mat)] <- 0
}

# --- B. TAX TABLE (FROM DEDICATED TAXONOMY FILE) ---
tax_file_dt <- fread(TAXONOMY, showProgress = FALSE)
setnames(tax_file_dt, tolower(names(tax_file_dt)))

tax_file_dt[, (target_id) := as.character(get(target_id))]
tax_file_dt <- unique(tax_file_dt, by = target_id)

if (anyDuplicated(tax_file_dt[[target_id]])) {
  warning(sprintf("⚠️ WARNING: Duplicate %s found in taxonomy file — keeping first occurrence only.", target_id))
}

tax_dt <- unique(tax_file_dt, by = target_id)

# Select available taxonomy columns
available_tax_cols <- base::intersect(tax_cols, names(tax_dt))
tax_mat <- as.matrix(tax_dt[, available_tax_cols, with = FALSE])
rownames(tax_mat) <- tax_dt[[target_id]]

# Re-align TAX table rows strictly with OTU matrix rows
tax_mat <- tax_mat[rownames(otu_mat), , drop = FALSE]

n_missing_tax <- sum(is.na(tax_mat[, 1]))
if (n_missing_tax > 0) {
  warning(sprintf("⚠️ WARNING: %d taxa present in OTU table have no taxonomy annotation.", n_missing_tax))
}

# --- C. SAMPLE DATA (STRICT RE-ALIGNMENT & SAFEGUARD) ---
sample_df <- as.data.frame(meta_dt)
rownames(sample_df) <- sample_df$sample_id
sample_df <- sample_df[colnames(otu_mat), , drop = FALSE]

if (anyNA(sample_df$sample_id)) {
  stop("❌ ERROR: Some samples present in the OTU matrix are missing from metadata after realignment.")
}

# ==========================================================================
# 4. Assemble and Save Object
# ==========================================================================
ps_final <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

saveRDS(ps_final, RDS)
message("✓ Phyloseq RDS object successfully created: ", RDS)