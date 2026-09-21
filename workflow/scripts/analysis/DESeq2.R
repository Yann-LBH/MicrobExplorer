# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : DESeq2.R
# PURPOSE : Differential Abundance & Expression Analysis with DESeq2
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
  library(readxl)
  library(arrow)
  #Libraries Bioconductor
  library(DESeq2)
})

source("workflow/scripts/utils/utils_io.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ==========================================================================
# Utility Function: Extract Contrast Results
# ==========================================================================
extract_results <- function(dds, contrast_vec, nom_contraste, extra_cols) {
  res <- results(dds, contrast = contrast_vec, cooksCutoff = FALSE)
  dt <- as.data.table(as.data.frame(res), keep.rownames = "Feature_ID")
  dt[, Comparison := nom_contraste]
  for (col in names(extra_cols)) dt[, (col) := extra_cols[[col]]]
  dt
}

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_DATA          <- as.character(snakemake@input[["data"]])
IN_METADATA      <- as.character(snakemake@input[["metadata"]])[1]

# Outputs
OUT_RDS           <- as.character(snakemake@output[["rds"]])[1]
OUT_PARQUET       <- as.character(snakemake@output[["parquet"]])[1]

# Controls and parameters
PARAM_CONTRAST      <- as.character(snakemake@params[["contrast"]])
PARAM_REF           <- as.character(snakemake@params[["ref"]])[1]
PARAM_SIZEFACTOR    <- as.character(snakemake@params[["sizefactor"]])[1] %||% "ratio"
PARAM_TEST          <- as.character(snakemake@params[["test"]])[1]       %||% "Wald"
PARAM_FITTYPE       <- as.character(snakemake@params[["fittype"]])[1]    %||% "parametric"

# Wildcards
WILDCARD_SOURCE        <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
if (is.null(IN_DATA) || length(IN_DATA) == 0 || !file.exists(IN_DATA[1])) {
  stop(sprintf("❌ Critical Error: Data input file '%s' does not exist.", IN_DATA[1]))
}

if (is.null(IN_METADATA) || !file.exists(IN_METADATA)) {
  stop(sprintf("❌ Critical Error: Metadata file '%s' does not exist.", IN_METADATA))
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading metadata...")
meta_dt <- load_metadata(IN_METADATA)
meta_dt[, date_real := as.Date(date, format = "%d/%m/%Y")]

if (!"group" %in% names(meta_dt)) {
  meta_dt[, group := name]
}

# Validation/Fallback pour la référence (PARAM_REF)
if (is.null(PARAM_REF) || PARAM_REF == "" || is.na(PARAM_REF) || !(PARAM_REF %in% meta_dt$group)) {
  if (!is.null(PARAM_REF) && PARAM_REF != "" && !is.na(PARAM_REF)) {
    warning("⚠️ WARNING: The 'ref' defined in config (", PARAM_REF, ") was not found in the Excel 'group' column. Falling back to default.\n")
  }
  PARAM_REF <- sort(meta_dt$group)[1]
}

message("INFO: Building count matrix...")
count_matrix <- build_deseq_count_matrix(IN_DATA, valid_ids = meta_dt$sample_id)
count_matrix <- count_matrix[rownames(count_matrix) != "KO_Unassigned", , drop = FALSE]

# Synchronisation des métadonnées (maintien sous forme de data.table)
meta_dt <- meta_dt[match(colnames(count_matrix), sample_id)]

if (nrow(meta_dt) == 0 || ncol(count_matrix) == 0) {
  stop("❌ Critical Error: Zero samples remaining after metadata alignment.")
}

message(sprintf("✓ Count matrix and metadata synchronized: %d samples, %d features.", 
                ncol(count_matrix), nrow(count_matrix)))

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PROCESSING FUNCTIONS
# ------------------------------------------------------------------------------

# Function 1: Group Condition Analysis (Group vs Ref & Pairwise Combos)
run_deseq_group_analyse <- function(count_matrix, meta_dt, PARAM_REF, PARAM_CONTRAST) {
  if (length(PARAM_CONTRAST) > 0 && !any(c("ref", "combo") %in% tolower(PARAM_CONTRAST))) return(NULL)

  col_data <- as.data.frame(meta_dt[, .(sample_id, group)])
  rownames(col_data) <- col_data$sample_id
  col_data <- col_data[colnames(count_matrix), , drop = FALSE]
  
  if (anyNA(col_data$sample_id)) {
    stop("❌ ERROR: Some samples in count_matrix are missing from metadata after realignment.")
  }
  col_data$group <- relevel(as.factor(col_data$group), ref = PARAM_REF)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ group)
  dds <- estimateSizeFactors(dds, type = PARAM_SIZEFACTOR)
  dds <- DESeq(dds, test = PARAM_TEST, fitType = PARAM_FITTYPE)

  groupes <- levels(col_data$group)

  dt_ref <- if (length(groupes) >= 2L && (length(PARAM_CONTRAST) == 0 || "ref" %in% PARAM_CONTRAST)) {
    rbindlist(lapply(groupes[groupes != PARAM_REF], function(g1) {
      extract_results(dds, c("group", g1, PARAM_REF), paste0(g1, "_vs_", PARAM_REF), list(Test_Group = g1, Ref_Group = PARAM_REF))
    }))
  } else NULL

  dt_combo <- if (length(groupes) >= 2L && (length(PARAM_CONTRAST) == 0 || "combo" %in% PARAM_CONTRAST)) {
    combos <- combn(groupes, 2, simplify = FALSE)
    rbindlist(lapply(combos, function(pair) {
      g2 <- pair[1]; g1 <- pair[2]
      extract_results(dds, c("group", g1, g2), paste0(g1, "_vs_", g2), list(Test_Group = g1, Ref_Group = g2))
    }))
  } else NULL

  message("✓ By group reference (ref ", PARAM_REF, ") : ", if (is.null(dt_ref)) 0 else nrow(dt_ref), " rows generated.")
  message("✓ By group combinations (all combos) : ", if (is.null(dt_combo)) 0 else nrow(dt_combo), " rows generated.")
  message("✓ Group analyses (REF + Combos) completed.")
  
  return(list(dds = dds, dt_ref = dt_ref, dt_combo = dt_combo))
}

# Function 2: Chronological Analysis (T vs T-1)
run_deseq_by_date <- function(count_matrix, meta_dt, PARAM_CONTRAST) {
  if (length(PARAM_CONTRAST) > 0 && !"date" %in% PARAM_CONTRAST) return(NULL)

  timeline <- unique(meta_dt[, .(date_real, date)])[order(date_real)]

  if (nrow(timeline) < 2L) {
    warning("Fewer than 2 distinct dates: skipping T vs T-1 contrast optimization.")
    return(invisible(NULL))
  }

  col_data <- as.data.frame(meta_dt[, .(sample_id, date, date_real, group)])
  rownames(col_data) <- col_data$sample_id
  col_data <- col_data[colnames(count_matrix), , drop = FALSE]
  
  if (anyNA(col_data$sample_id)) {
    stop("❌ ERROR: Some samples in count_matrix are missing from metadata after realignment.")
  }
  col_data$date_group <- as.factor(col_data$date)
  col_data$group      <- as.factor(col_data$group)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ group + date_group)
  dds <- estimateSizeFactors(dds, type = PARAM_SIZEFACTOR)
  dds <- DESeq(dds, test = PARAM_TEST, fitType = PARAM_FITTYPE)

  dt_date <- rbindlist(lapply(
    2:nrow(timeline),
    function(i) {
      t_curr <- as.character(timeline$date[i])
      t_prev <- as.character(timeline$date[i - 1L])
      extract_results(
        dds,
        c("date_group", t_curr, t_prev),
        paste0(t_curr, "_vs_", t_prev),
        list(Date_Test = t_curr, Date_Ref = t_prev)
      )
    }
  ))

  message("✓ By timeline (T vs T-1) : ", nrow(dt_date), " rows generated.")
  return(list(dt = dt_date, dds = dds))
}

# ------------------------------------------------------------------------------
# 6. EXECUTION CORE & GRAPHICS GENERATION
# ------------------------------------------------------------------------------
res_group <- run_deseq_group_analyse(count_matrix, meta_dt, PARAM_REF, PARAM_CONTRAST) 
res_date  <- run_deseq_by_date(count_matrix, meta_dt, PARAM_CONTRAST)

models_list <- list()
dt_list     <- list()

if (!is.null(res_group) && !is.null(res_group$dt_ref) && nrow(res_group$dt_ref) > 0) {
  models_list$ref <- res_group$dds
  dt_list$ref     <- res_group$dt_ref[, Contrast_Type := "ref"]
}

if (!is.null(res_group) && !is.null(res_group$dt_combo) && nrow(res_group$dt_combo) > 0) {
  models_list$combo <- res_group$dds
  dt_list$combo     <- res_group$dt_combo[, Contrast_Type := "combo"]
}

if (!is.null(res_date) && !is.null(res_date$dt) && nrow(res_date$dt) > 0) {
  models_list$date <- res_date$dds
  dt_list$date     <- res_date$dt[, Contrast_Type := "date"]
}

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (length(dt_list) > 0) {
  # Compile and save Master OUT_RDS
  master_rds <- list(
    models = models_list,
    dt     = dt_list
  )
  saveRDS(master_rds, OUT_RDS)

  # Compile and export Parquet
  master_parquet <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  arrow::write_parquet(master_parquet, OUT_PARQUET)

  message("✓ Success exports written:")
  message("  - RDS     : ", OUT_RDS)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  # Save empty objects if no contrast was generated
  saveRDS(list(models = list(), dt = list()), OUT_RDS)
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: No contrast results generated. Empty Parquet written.")
}