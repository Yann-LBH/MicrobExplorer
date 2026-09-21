# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Phyloseq.R
# PURPOSE : Build Phyloseq RDS object from Abundance Data, Metadata & Taxonomy
# AUTHOR  : Yann Le Bihan
# DATE    : 2025-12-01
# LINK    : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. METADATA & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Libraries CRAN
  library(data.table)
  library(readxl)
  # Libraries Bioconductor
  library(phyloseq)
})

source("workflow/scripts/utils/utils_io.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_DATA             <- as.character(snakemake@input[["data"]])
IN_METADATA         <- as.character(snakemake@input[["metadata"]])[1]
IN_TAXONOMY         <- as.character(snakemake@input[["taxonomy"]])[1]

# Outputs
OUT_RDS              <- as.character(snakemake@output[["rds"]])[1]

# Parameters
PARAM_STAND_COL  <- tolower(as.character(snakemake@params[["stand_col"]]))[1]
PARAM_RANK_KEGG  <- as.character(snakemake@params[["rank"]])[1]

# Wildcards
WILDCARD_SOURCE  <- if (!is.null(snakemake@params[["source"]]) && !is.na(snakemake@params[["source"]])) {
  tolower(as.character(snakemake@params[["source"]]))[1]
} else if (!is.null(snakemake@wildcards[["source"]]) && !is.na(snakemake@wildcards[["source"]])) {
  tolower(as.character(snakemake@wildcards[["source"]]))[1]
} else {
  NA_character_
}

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
if (is.null(IN_DATA) || length(IN_DATA) == 0 || !file.exists(IN_DATA[1])) {
  stop(sprintf("❌ Critical Error: Data input path '%s' does not exist.", IN_DATA[1]))
}

if (is.null(IN_METADATA) || !file.exists(IN_METADATA)) {
  stop(sprintf("❌ Critical Error: Metadata file '%s' does not exist.", IN_METADATA))
}

if (is.null(IN_TAXONOMY) || !file.exists(IN_TAXONOMY)) {
  stop(sprintf("❌ Critical Error: Taxonomy file '%s' does not exist.", IN_TAXONOMY))
}

if (is.null(PARAM_STAND_COL) || is.na(PARAM_STAND_COL) || PARAM_STAND_COL == "") {
  stop("❌ Critical Error: 'stand_col' parameter is missing or empty in Snakemake config.")
}

if (is.na(WILDCARD_SOURCE) || !(WILDCARD_SOURCE %in% c("kegg", "reads", "contigs"))) {
  stop(sprintf(
    "❌ Critical Error: Invalid or missing SOURCE ['%s']. Expected 'kegg', 'reads', or 'contigs'. Check params/wildcards in Snakefile.",
    WILDCARD_SOURCE
  ))
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
# Dynamic Configuration Based on the Source
is_kegg  <- grepl("kegg", WILDCARD_SOURCE, ignore.case = TRUE)
use_rank <- is_kegg && !is.null(PARAM_RANK_KEGG) && !is.na(PARAM_RANK_KEGG) && !(PARAM_RANK_KEGG %in% c("", "NA", "all"))
RANK     <- if (use_rank) PARAM_RANK_KEGG else NULL

if (use_rank) {
  message(sprintf("ℹ️ [KEGG Mode] Specific rank parameter enabled: '%s'", RANK))
} else if (is_kegg) {
  message("ℹ️ [KEGG Mode] No specific rank provided. Using all KEGG ranks.")
} else {
  message(sprintf("ℹ️ [Mode %s] 'rank' parameter ignored (only active for KEGG).", toupper(WILDCARD_SOURCE)))
}

if (is_kegg) {
  target_id <- "kegg_id"
  tax_cols  <- c("ec_number", "level_1", "level_2", "level_3", "gene_description")
  
  if (use_rank && !(RANK %in% tax_cols) && RANK != "kegg_id") {
    warning(sprintf("⚠️ WARNING: The rank '%s' is not in the valid KEGG columns (%s).", 
                    RANK, paste(tax_cols, collapse = ", ")))
  }
} else if (grepl("reads", WILDCARD_SOURCE, ignore.case = TRUE)) {
  target_id <- "read_id" 
  tax_cols  <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
} else if (grepl("contigs", WILDCARD_SOURCE, ignore.case = TRUE)) {
  target_id <- "contig_id"
  tax_cols  <- c("domain", "phylum", "class", "order", "family", "genus", "species")
} else {
  stop(sprintf("❌ Critical Error: Unknown SOURCE value [%s]. Expected 'kegg', 'reads', or 'contigs'.", WILDCARD_SOURCE))
}

message("INFO: Loading metadata & abundance data...")
meta_dt <- load_metadata(IN_METADATA)

all_data <- load_tsv_dir_dynamic(
  paths       = IN_DATA, 
  meta_dt     = meta_dt, 
  select_cols = c(target_id, PARAM_STAND_COL)
)

all_data[[target_id]] <- as.character(all_data[[target_id]])

# Cleaning Up Missing or Empty Values
valid_mask <- !is.na(all_data[[target_id]]) & 
              all_data[[target_id]] != "" & 
              all_data[[target_id]] != "NA"
all_data   <- all_data[valid_mask]
all_data[[PARAM_STAND_COL]] <- as.numeric(all_data[[PARAM_STAND_COL]])

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PHYLOSEQ COMPONENTS
# ------------------------------------------------------------------------------
message("INFO: Loading taxonomy reference table...")
taxo_ref <- load_taxonomy_table(
  taxonomy_path = IN_TAXONOMY, 
  source        = WILDCARD_SOURCE
)

tax_dt   <- data.table::copy(taxo_ref$dt)
join_col <- taxo_ref$join_col
tax_dt   <- tax_dt[!is.na(tax_dt[[join_col]]) & tax_dt[[join_col]] != ""]

pathway_levels <- c("level_1", "level_2", "level_3")

is_pathway_aggregation <- use_rank && RANK %in% pathway_levels

if (is_pathway_aggregation) {
  # --- CASE A: AGGREGATION BY PATHWAY ---
  message(sprintf("🔄 [KEGG] Aggregating counts at rank: '%s'", RANK))
  
  # 1. Join counts data with taxonomy annotations
  all_data_merged <- merge(all_data, tax_dt, by.x = target_id, by.y = join_col, allow.cartesian = TRUE)
  all_data_merged <- all_data_merged[!is.na(get(RANK)) & get(RANK) != "" & get(RANK) != "Unclassified"]
  
  # 2. Build aggregated OTU matrix by target pathway rank
  formula_str <- as.formula(paste(RANK, "~ sample_id"))
  otu_dt <- dcast(all_data_merged, formula_str, value.var = PARAM_STAND_COL, fun.aggregate = sum, fill = 0)
  otu_mat <- as.matrix(otu_dt, rownames = RANK)
  
  # 3. Identify parent levels to keep up to current RANK
  idx_rank      <- which(pathway_levels == RANK)
  ranks_to_keep <- pathway_levels[1:idx_rank]
  
  # 4. Aggregated taxonomy table grouped uniquely by RANK
  # Parent ranks and KOs are concatenated if a pathway appears under multiple parent categories
  tax_dt_rank <- tax_dt[!is.na(get(RANK)) & get(RANK) != "", lapply(.SD, function(x) {
    u_vals <- unique(na.omit(x[x != ""]))
    if (length(u_vals) == 0) return("Unassigned")
    paste(sort(u_vals), collapse = "; ")
  }), by = RANK, .SDcols = c(setdiff(ranks_to_keep, RANK), join_col)]
  
  # Rename the concatenated join column to 'kegg_ids'
  setnames(tax_dt_rank, old = join_col, new = "kegg_ids")
  
  # Reorder columns to follow hierarchy order
  cols_order <- c(ranks_to_keep, "kegg_ids")
  tax_mat_raw <- as.matrix(tax_dt_rank[, cols_order, with = FALSE])
  rownames(tax_mat_raw) <- tax_dt_rank[[RANK]]
} else {
  # --- CASE B: GENE LEVEL / KO / READ / CONTIG ---
  formula_str <- as.formula(paste(target_id, "~ sample_id"))
  otu_dt      <- dcast(all_data, formula_str, value.var = PARAM_STAND_COL, fun.aggregate = sum, fill = 0)
  otu_mat     <- as.matrix(otu_dt, rownames = target_id)

  if (is_kegg) {
    ranks_to_use  <- intersect(tax_cols, colnames(tax_dt))
    tax_dt_unique <- tax_dt[, lapply(.SD, function(x) {
      u_vals <- unique(na.omit(x[x != ""]))
      if (length(u_vals) == 0) return("Unassigned")
      paste(u_vals, collapse = "; ")
    }), by = join_col, .SDcols = ranks_to_use]
    
    tax_mat_raw <- as.matrix(tax_dt_unique[, ranks_to_use, with = FALSE])
    rownames(tax_mat_raw) <- tax_dt_unique[[join_col]]
  } else {
    ranks_to_use <- intersect(tax_cols, colnames(tax_dt))
    tax_dt_unique <- unique(tax_dt, by = join_col)
    tax_mat_raw  <- as.matrix(tax_dt_unique[, ranks_to_use, with = FALSE])
    rownames(tax_mat_raw) <- tax_dt_unique[[join_col]]
  }
}

# --- ALIGNMENT AND SECURING OF MATRICES ---
mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0

feature_ids <- rownames(otu_mat)
matched_idx <- match(feature_ids, rownames(tax_mat_raw))

tax_mat <- tax_mat_raw[matched_idx, , drop = FALSE]
rownames(tax_mat) <- feature_ids

n_missing_tax <- sum(is.na(tax_mat[, 1]))
if (n_missing_tax > 0) {
  warning(sprintf("⚠️ WARNING: %d elements in the OTU table do not have taxonomy annotations in the taxonomy file.", n_missing_tax))
}

tax_mat[is.na(tax_mat)] <- "Unclassified"

# --- SAMPLES (METADATA ALIGNMENT) ---
sample_df           <- as.data.frame(meta_dt)
rownames(sample_df) <- sample_df$sample_id
sample_df           <- sample_df[colnames(otu_mat), , drop = FALSE]

if (anyNA(sample_df$sample_id)) {
  stop("❌ Critical Error: Some samples present in the OTU matrix are missing from metadata after realignment.")
}

# ------------------------------------------------------------------------------
# 6. GRAPHICS GENERATION & OBJECT ASSEMBLY
# ------------------------------------------------------------------------------
message("INFO: Assembling Phyloseq object...")
ps_final <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (!is.null(ps_final) && nsamples(ps_final) > 0) {
  saveRDS(ps_final, OUT_RDS)

  message("✓ Success exports written:")
  message("  - RDS     : ", OUT_RDS)
} else {
  saveRDS(ps_final, OUT_RDS)
  warning("⚠️ WARNING: No valid effects evaluated. Empty RDS file created.")
}