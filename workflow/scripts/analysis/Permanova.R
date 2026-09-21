# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Permanova.R
# PURPOSE : PERMANOVA (adonis2) and Multivariate Dispersion (betadisper)
# AUTHOR  : Yann Le Bihan
# DATE    : 2026-09-03
# LINK    : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. METADATA & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Libraries CRAN
  library(data.table)
  library(rlang)
  library(arrow)
  library(compositions) #CLR
  library(vegan)  # PERMANOVA, betadisper, decostand

})

source("workflow/scripts/utils/utils_io.R")
source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_DATA            <- as.character(snakemake@input[["data"]])
IN_METADATA        <- as.character(snakemake@input[["metadata"]])[1]

# Outputs
OUT_PDF             <- as.character(snakemake@output[["pdf"]])[1]
OUT_XLSX            <- as.character(snakemake@output[["xlsx"]])[1]
OUT_PARQUET         <- as.character(snakemake@output[["parquet"]])[1]

# Shared plot features
PARAM_SHARED          <- snakemake@params[["shared"]]
PARAM_PDF_SIZE        <- as.numeric(PARAM_SHARED$pdf_size) %||% c(14, 12)

# Parameters
PARAM_STAND_COL       <- tolower(as.character(snakemake@params[["stand_col"]]))[1]
PARAM_EFFECT          <- as.character(snakemake@params[["effect"]])
PARAM_DISTANCE_METHOD <- tolower(as.character(snakemake@params[["distance_method"]]))[1]
PARAM_PERMUTATIONS    <- as.numeric(snakemake@params[["permutation"]])[1]
PARAM_USE_CLR         <- isTRUE(snakemake@params[["use_clr"]])

# Wildcards
WILDCARD_SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
if (is.null(IN_DATA) || length(IN_DATA) == 0 || !file.exists(IN_DATA[1])) {
  stop(sprintf("❌ Critical Error: Data input path '%s' does not exist.", IN_DATA[1]))
}

if (is.null(IN_METADATA) || !file.exists(IN_METADATA)) {
  stop(sprintf("❌ Critical Error: Metadata file '%s' does not exist.", IN_METADATA))
}

if (is.null(PARAM_STAND_COL) || is.na(PARAM_STAND_COL) || PARAM_STAND_COL == "") {
  stop("❌ Critical Error: 'stand_col' parameter is missing or empty in Snakemake config.")
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading metadata...")
meta_dt <- load_metadata(IN_METADATA)

message("INFO: Loading abundance TSV files...")
counts_dt <- load_tsv_dir_dynamic(
  paths   = IN_DATA,
  meta_dt = meta_dt
)

# Explicit identification of the 'feature_id' column
feature_col <- intersect(c("read_id", "contig_id", "kegg_id"), colnames(counts_dt))[1]
if (is.na(feature_col)) {
  stop(sprintf(
    "❌ Critical Error: No known feature ID column found. Available columns: %s",
    paste(colnames(counts_dt), collapse = ", ")
  ))
}

# Checking the value column (stand_col)
value_col <- PARAM_STAND_COL
if (!value_col %in% colnames(counts_dt)) {
  stop(sprintf(
    "❌ Critical Error: Column '%s' not found in loaded data. Available columns: %s",
    value_col, paste(colnames(counts_dt), collapse = ", ")
  ))
}

# STRICT INTEGRITY CHECK: Aggregation Verification
dup_rows <- duplicated(counts_dt[, c("sample_id", feature_col), with = FALSE])

if (any(dup_rows)) {
  dup_check            <- counts_dt[dup_rows, ]
  all_affected_samples <- unique(dup_check$sample_id)
  samples_str          <- paste(all_affected_samples, collapse = ", ")
  
  stop(sprintf(
    "❌ Critical Error: Non-aggregated data detected! Multiple entries found for the same [sample_id, %s] pair.\n Total duplicate entries: %d\n All affected samples (%d total): %s\n",
    feature_col, 
    nrow(dup_check), 
    length(all_affected_samples), 
    samples_str
  ))
}

# Switch from LONG format to WIDE format (Samples x Features)
dcast_formula <- as.formula(sprintf("sample_id ~ %s", feature_col))
matrix_dt <- dcast(
  counts_dt, 
  dcast_formula, 
  value.var = value_col, 
  fill = 0
)

# Convert to a numeric matrix with `sample_id` as row names
ready_data           <- as.matrix(matrix_dt[, -1, with = FALSE])
rownames(ready_data) <- matrix_dt$sample_id

# Strict alignment of metadata
metadata           <- as.data.frame(meta_dt)
rownames(metadata) <- metadata$sample_id
common_samples     <- intersect(rownames(ready_data), rownames(metadata))

ready_data <- ready_data[common_samples, , drop = FALSE]
metadata   <- metadata[common_samples, , drop = FALSE]

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PROCESSING FUNCTIONS
# ------------------------------------------------------------------------------
if (PARAM_USE_CLR) {
  message("INFO: Applying CLR transformation (with pseudocount +1)...")
  ready_data      <- as.matrix(clr(ready_data + 1))
  PARAM_DISTANCE_METHOD <- "euclidean"
} else {
  message("INFO: Applying relative abundance normalization (TSS)...")
  ready_data      <- vegan::decostand(ready_data, method = "total")
}

# ------------------------------------------------------------------------------
# 6. EXECUTION CORE & GRAPHICS GENERATION
# ------------------------------------------------------------------------------
message("INFO: Running PERMANOVA and Dispersion analysis...")

dist_matrix <- vegan::vegdist(ready_data, method = PARAM_DISTANCE_METHOD)

results_list <- list()

with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  page_count <- 0

  for (eff in PARAM_EFFECT) {
    if (!eff %in% colnames(metadata)) {
      warning(sprintf("⚠️ WARNING: The '%s' effect is not present in metadata. Ignored.", eff))
      next
    }

    # Handling NA values in the variable of interest
    valid_idx <- !is.na(metadata[[eff]])
    if (sum(valid_idx) < 3) {
      warning(sprintf("⚠️ WARNING: Not enough valid samples (< 3) for effect '%s'. Ignored.", eff))
      next
    }

    # Dynamic subsetting if NA values exist (otherwise, retain 100% of the data)
    sub_meta <- metadata[valid_idx, , drop = FALSE]
    
    if (all(valid_idx)) {
      sub_dist <- dist_matrix
    } else {
      message(sprintf("INFO: Removing %d sample(s) with NA for effect '%s'", sum(!valid_idx), eff))
      # Extraction du sous-ensemble de la matrice de distance
      sub_dist <- as.dist(as.matrix(dist_matrix)[valid_idx, valid_idx])
    }

    message(sprintf("Processing effect: %s", eff))
    
    # PERMANOVA (adonis2)
    formula_permanova <- as.formula(sprintf("sub_dist ~ %s", eff))
    perm_res          <- vegan::adonis2(
      formula_permanova, 
      data         = sub_meta,
      permutations = PARAM_PERMUTATIONS
    )
    
    df_permanova <- as.data.frame(perm_res)
    df_permanova <- cbind(
      Test     = paste("PERMANOVA -", eff),
      Variable = rownames(df_permanova),
      df_permanova
    )

    # Dispersion Homogeneity (betadisper)
    disper      <- vegan::betadisper(sub_dist, sub_meta[[eff]])

    plot(disper, main = sprintf("Multivariate Dispersion: %s (%s)", eff, WILDCARD_SOURCE))
    render_page(recordPlot())

    disper_anova  <- as.data.frame(anova(disper))
    df_dispersion <- cbind(
      Test     = paste("DISPERSION -", eff),
      Variable = rownames(disper_anova),
      disper_anova
    )

    # Save the results
    results_list[[paste0("permanova_", eff)]]  <- df_permanova
    results_list[[paste0("dispersion_", eff)]] <- df_dispersion
    page_count <- page_count + 1
  }
   if (page_count == 0) {
    render_fallback("No valid effects evaluated for PERMANOVA.")
  }
})
# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (length(results_list) > 0) {
  export_final <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

  # Cleaning up column names for Parquet compatibility (e.g., Pr(>F) -> Pr_F)
  setnames(export_final, old = names(export_final), new = gsub("[^A-Za-z0-9_]", "_", names(export_final)))

  writexl::write_xlsx(export_final, OUT_XLSX)
  arrow::write_parquet(export_final, OUT_PARQUET)

  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Excel     :", OUT_XLSX)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  writexl::write_xlsx(data.frame(), OUT_XLSX)
  arrow::write_parquet(data.frame(), OUT_PARQUET)
  warning("⚠️ WARNING: No valid effects evaluated. Empty Parquet file created.")
}