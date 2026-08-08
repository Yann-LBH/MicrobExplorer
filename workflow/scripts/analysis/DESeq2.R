################################################################################
# Project : "MicrobExplorer"
# Script  : "Analysis : DESeq2 on raw counts"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

suppressPackageStartupMessages({
  # Libraries CRAN
  library(data.table)
  library(readxl)
  #Libraries Bioconductor
  library(DESeq2)
  library(arrow)
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
PARQUET  <- as.character(snakemake@output[["parquet"]])[1]

# Controls and parameters
CONTRAST_LIST <- tolower(as.character(snakemake@params[["contrast"]]))
REF       <- as.character(snakemake@params[["ref"]])[1]
SIZEFACTOR <- as.character(snakemake@params[["sizefactor"]])[1]
TEST      <- as.character(snakemake@params[["test"]])[1]
FITTYPE   <- as.character(snakemake@params[["fittype"]])[1]

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ==========================================================================
# 1. Loading Metadata and Files
# ==========================================================================

# 1. Metadata loading
meta_dt <- load_metadata(METADATA)
meta_dt[, Date_Real := as.Date(date, format = "%d/%m/%Y")]

if (!"group" %in% names(meta_dt)) {
  meta_dt[, group := name]
}

# Validate if the REF from config exists in your Excel "group" column
if (is.null(REF) || REF == "" || is.na(REF) || !(REF %in% meta_dt$group)) {
  if (!is.null(REF) && REF != "" && !is.na(REF)) {
    warning("⚠️ WARNING: The 'ref' defined in config (", REF, ") was not found in the Excel 'group' column. Falling back to default.\n")
  }
  REF <- sort(meta_dt$group)[1]
}

# 2. Chargement du count_matrix (sécurités déportées dans utils_io.R)
count_matrix <- build_deseq_count_matrix(DATA, valid_ids = meta_dt$sample_id)

# 3. Synchronize metadata (Keep it as data.table to prevent downstream "." query crashes)
meta_dt <- meta_dt[match(colnames(count_matrix), sample_id)]

# ==========================================================================
# TAXONOMY loading and merging
# ==========================================================================
load_taxonomy_table <- function(taxonomy_path, source) {
  if (!file.exists(taxonomy_path)) {
    stop(sprintf("❌ ERROR: Taxonomy file not found: %s", taxonomy_path))
  }
  dt_taxo <- fread(taxonomy_path, showProgress = FALSE)
  setnames(dt_taxo, tolower(names(dt_taxo)))

  if (source == "reads") {
    join_col <- "tax_id"
    tax_ranks <- c("tax_id", "scientific_name", "domain", "kingdom",
                   "phylum", "class", "order", "family", "genus", "species")
  } else if (source == "contigs") {
    join_col <- "contig_id"
    tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  } else {
    join_col <- "ko"
    tax_ranks <- c("ec_number", "level_1", "level_2", "level_3", "gene_description")
  }

  if (!join_col %in% names(dt_taxo)) {
    stop(sprintf(
      "❌ ERROR: Expected join column '%s' not found in taxonomy file. Available columns: %s",
      join_col, paste(names(dt_taxo), collapse = ", ")
    ))
  }

  dt_taxo[, (join_col) := as.character(get(join_col))]

  if (anyDuplicated(dt_taxo[[join_col]])) {
    n_dup <- sum(duplicated(dt_taxo[[join_col]]))
    warning(sprintf(
      "⚠️ WARNING: %d duplicated '%s' found in taxonomy file — keeping first occurrence only.",
      n_dup, join_col
    ))
  }
  dt_taxo <- unique(dt_taxo, by = join_col)

  tax_ranks_present <- intersect(tax_ranks, names(dt_taxo))
  list(dt = dt_taxo[, c(join_col, tax_ranks_present), with = FALSE],
       join_col = join_col,
       tax_ranks = tax_ranks_present)
}

annotate_with_taxonomy <- function(dt_results, taxo_ref, id_col_results = "Feature_ID") {
  dt_results[, (id_col_results) := as.character(get(id_col_results))]

  dt_annotated <- merge(
    dt_results, taxo_ref$dt,
    by.x = id_col_results, by.y = taxo_ref$join_col, all.x = TRUE
  )

  for (col in taxo_ref$tax_ranks) {
    n_na <- sum(is.na(dt_annotated[[col]]))
    if (n_na > 0) {
      set(dt_annotated, i = which(is.na(dt_annotated[[col]])), j = col, value = "Unclassified")
    }
  }
  dt_annotated
}

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

# ==========================================================================
# 2. Condition Analysis (Group vs Reference) & (Combo)
# ==========================================================================
run_deseq_group_analyse <- function(count_matrix, meta_dt, REF, CONTRAST_LIST) {
  if (length(CONTRAST_LIST) > 0 && !any(c("ref", "combo") %in% CONTRAST_LIST)) return(NULL)

  col_data <- as.data.frame(meta_dt[, .(sample_id, group)])
  rownames(col_data) <- col_data$sample_id
  col_data <- col_data[colnames(count_matrix), , drop = FALSE]
  if (anyNA(col_data$sample_id)) {
    stop("❌ ERROR: Some samples in count_matrix are missing from metadata after realignment.")
  }
  col_data$group <- relevel(as.factor(col_data$group), ref = REF)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ group)
  dds <- estimateSizeFactors(dds, type = SIZEFACTOR)
  dds <- DESeq(dds, test = TEST, fitType = FITTYPE)

  groupes <- levels(col_data$group)

  dt_ref <- if (length(groupes) >= 2L && (length(CONTRAST_LIST) == 0 || "ref" %in% CONTRAST_LIST)) {
    rbindlist(lapply(groupes[groupes != REF], function(g1) {
      extract_results(dds, c("group", g1, REF), paste0(g1, "_vs_", REF), list(Test_Group = g1, Ref_Group = REF))
    }))
  } else NULL

  # Pairwise combos (executed if requested or by default, and if >= 2 groups)
  dt_combo <- if (length(groupes) >= 2L && (length(CONTRAST_LIST) == 0 || "combo" %in% CONTRAST_LIST)) {
    combos <- combn(groupes, 2, simplify = FALSE)
    rbindlist(lapply(combos, function(pair) {
      g2 <- pair[1]; g1 <- pair[2]
      extract_results(dds, c("group", g1, g2), paste0(g1, "_vs_", g2), list(Test_Group = g1, Ref_Group = g2))
    }))
  } else NULL
  
  message("✓ By group reference (ref ", REF, ") : ", if (is.null(dt_ref)) 0 else nrow(dt_ref), " rows generated.")
  message("✓ By group combinations (all combos) : ", if (is.null(dt_combo)) 0 else nrow(dt_combo), " rows generated.")
  message("✓ Group analyses (REF + Combos) completed.")
  return(list(dds = dds, dt_ref = dt_ref, dt_combo = dt_combo))
}

# ==========================================================================
# 3. Chronological Analysis (T vs T-1)
# ==========================================================================
run_deseq_by_date <- function(count_matrix, meta_dt, CONTRAST_LIST) {
  if (length(CONTRAST_LIST) > 0 && !"date" %in% CONTRAST_LIST) return(NULL)

  timeline <- unique(meta_dt[, .(Date_Real, date)])[order(Date_Real)]

  if (nrow(timeline) < 2L) {
    warning("Fewer than 2 distinct dates: skipping T vs T-1 contrast optimization.")
    return(invisible(NULL))
  }

  col_data <- as.data.frame(meta_dt[, .(sample_id, date, Date_Real)])
  rownames(col_data) <- col_data$sample_id
  col_data <- col_data[colnames(count_matrix), , drop = FALSE]
  if (anyNA(col_data$sample_id)) {
    stop("❌ ERROR: Some samples in count_matrix are missing from metadata after realignment.")
  }
  col_data$Date_Group <- as.factor(col_data$date)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~Date_Group)
  dds <- estimateSizeFactors(dds, type = SIZEFACTOR)
  dds <- DESeq(dds, test = TEST, fitType = FITTYPE)

  dt_date <- rbindlist(lapply(
    2:nrow(timeline),
    function(i) {
      t_curr <- as.character(timeline$date[i])
      t_prev <- as.character(timeline$date[i - 1L])
      extract_results(
        dds,
        c("Date_Group", t_curr, t_prev),
        paste0(t_curr, "_vs_", t_prev),
        list(Date_Test = t_curr, Date_Ref = t_prev)
      )
    }
  ))

  message("✓ By timeline (T vs T-1) : ", nrow(dt_date), " rows generated.")
  return(list(dt = dt_date, dds = dds))
}

# ==========================================================================
# Execution Core
# ==========================================================================

# 1. Run all analyses and return both the results table AND the dds object from each
res_group    <- run_deseq_group_analyse(count_matrix, meta_dt, REF, CONTRAST_LIST) 
res_date   <- run_deseq_by_date(count_matrix, meta_dt, CONTRAST_LIST)

models_list  <- list()
dt_list      <- list()

# Dynamic detection of available results and appending them to the master lists
if (!is.null(res_group) && !is.null(res_group$dt_ref)) {
  models_list$ref <- res_group$dds
  dt_list$ref     <- res_group$dt_ref[, Contrast_Type := "ref"]
}

if (!is.null(res_group) && !is.null(res_group$dt_combo)) {
  models_list$combo <- res_group$dds
  dt_list$combo     <- res_group$dt_combo[, Contrast_Type := "combo"]
}

if (!is.null(res_date) && !is.null(res_date$dt)) {
  models_list$date <- res_date$dds
  dt_list$date     <- res_date$dt[, Contrast_Type := "date"]
}

# 2. Annotate all result tables in dt_list with taxonomy BEFORE saving
if (length(dt_list) > 0) {
  taxo_ref <- load_taxonomy_table(TAXONOMY, SOURCE)
  dt_list <- lapply(dt_list, annotate_with_taxonomy, taxo_ref = taxo_ref)
}

# 3. Compile EVERYTHING into the Master RDS file
master_rds <- list(
  models = models_list,
  dt     = dt_list
)
saveRDS(master_rds, RDS)

# 4. Compile and bind for global Parquet output
if (length(dt_list) > 0) {
  master_parquet <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  write_parquet(master_parquet, PARQUET)
} else {
  warning("⚠️ No contrast results generated to write to Parquet.")
}

message("✅ Master RDS (Models + Stats) and Parquet files compiled successfully.")