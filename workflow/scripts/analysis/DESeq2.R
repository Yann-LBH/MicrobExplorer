################################################################################
# Project : "MicrobExplorer"
# Script  : "Analysis : DESeq2"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(readxl)
library(DESeq2)
library(arrow)

DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])
RDS      <- as.character(snakemake@output[["rds"]])
PARQUET  <- as.character(snakemake@output[["parquet"]])

# Controls and parameters
CONTRASTS <- tolower(as.character(snakemake@params[["contrasts"]]))[1]
REF       <- as.character(snakemake@params[["ref"]])[1] # Case-sensitive matching (e.g., "TD1")

# ==========================================================================
# 1. Loading Metadata and Files
# ==========================================================================
# All script comments are provided in English as requested.

meta_dt <- as.data.table(read_xlsx(METADATA))
meta_dt[, Date_Real := as.Date(date, format = "%d/%m/%Y")]

# ✅ FIXED: Use your new Excel column "group" directly. Fallback to "name" only if missing.
if (!"group" %in% names(meta_dt)) {
  meta_dt[, group := name]
}

# Validate if the REF from config exists in your Excel "group" column
if (is.null(REF) || REF == "" || is.na(REF) || !(REF %in% meta_dt$group)) {
  if (!is.null(REF) && REF != "" && !is.na(REF)) {
    cat("⚠️ WARNING: The 'ref' defined in config (", REF, ") was not found in the Excel 'group' column. Falling back to default.\n")
  }
  REF <- sort(meta_dt$group)[1]
}

# Cross-load files matching existing sample ids
raw_list <- lapply(DATA, function(f) {
  file_name <- basename(f)
  
  matched_sample <- meta_dt$sample_id[sapply(meta_dt$sample_id, function(sid) grepl(sid, file_name))]

  if (length(matched_sample) == 0 || is.na(matched_sample[1]) || matched_sample[1] == "") {
    stop(paste0(
      "\n❌ ERROR: No matching sample_id found in metadata for file: '", file_name, "'\n",
      "Please check if this sample is declared in your metadata.xlsx or verify the filename."
    ))
  }
  
  dt <- fread(f, showProgress = FALSE)
  if (nrow(dt) == 0) return(NULL)

  setnames(dt, tolower(names(dt)))
  
  # Identify the ID column
  current_id <- base::intersect(c("read_id", "contig_id", "ko", "kegg"), names(dt))[1]
  if (is.na(current_id)) return(NULL)
  
  # Identify the count column and enforce numeric representation
  abundance_col <- base::intersect(c("count", "read_mapped"), names(dt))[1]
  if (is.na(abundance_col)) {
    abundance_col <- names(dt)[ncol(dt)]
  }
  dt[, (abundance_col) := lapply(.SD, as.numeric), .SDcols = abundance_col]
  
  setnames(dt, abundance_col, matched_sample[1])
  dt <- dt[, c(current_id, matched_sample[1]), with = FALSE]
  
  list(dt = dt, sample_id = matched_sample[1])
})

raw_list <- Filter(Negate(is.null), raw_list)

if (length(raw_list) == 0) {
  stop("🚨 Step error: raw_list is empty. No valid sample data tables were loaded for DESeq2 analysis.")
}

possible_ids <- base::intersect(c("contig_id", "read_id", "ko", "kegg"), names(raw_list[[1]]$dt))[1]
if (is.na(possible_ids)) {
  possible_ids <- "read_id"
}

# Merge count matrices
count_data <- Reduce(
  function(a, b) merge(a, b, by = possible_ids, all = FALSE),
  lapply(raw_list, `[[`, "dt")
)

count_matrix <- as.matrix(count_data[, !possible_ids, with = FALSE])
rownames(count_matrix) <- count_data[[possible_ids]]
count_matrix <- round(count_matrix)

# Synchronize metadata (Keep it as data.table to prevent downstream "." query crashes)
meta_dt <- meta_dt[sample_id %in% colnames(count_matrix)]

rm(raw_list)

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
# 2. Condition Analysis (Group vs Reference)
# ==========================================================================
run_deseq_by_name_ref <- function(count_matrix, meta_dt, REF, RDS, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, group)])
  rownames(col_data) <- col_data$sample_id
  
  col_data$group <- relevel(as.factor(col_data$group), ref = REF)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")

  groupes <- levels(col_data$group)

  results_dt <- rbindlist(lapply(
    groupes[groupes != REF],
    function(g1) {
      extract_results(
        dds,
        c("group", g1, REF),
        paste0(g1, "_vs_", REF),
        list(Test_Group = g1, Ref_Group = REF)
      )
    }
  ))
  
  return(list(dt = results_dt, dds = dds))
  return(results_dt)
  message("✓ By group reference (ref ", REF, ") : ", nrow(results_dt), " rows generated.")
}

# ==========================================================================
# 3. Chronological Analysis (T vs T-1)
# ==========================================================================
run_deseq_by_date <- function(count_matrix, meta_dt, RDS, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, date, Date_Real)])
  rownames(col_data) <- col_data$sample_id
  col_data$Date_Group <- as.factor(col_data$date)

  # ✅ INDEPENDENT: Continues to use Date_Group for timeline analysis
  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~Date_Group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")

  timeline <- unique(meta_dt[, .(Date_Real, date)])[order(Date_Real)]

  if (nrow(timeline) < 2L) {
    warning("Fewer than 2 distinct dates: skipping T vs T-1 contrast optimization.")
    return(invisible(NULL))
  }

  results_dt <- rbindlist(lapply(
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

  return(list(dt = results_dt, dds = dds))
  return(results_dt)
  message("✓ By timeline (T vs T-1) : ", nrow(results_dt), " rows generated.")
}

# ==========================================================================
# 4. Pairwise Combination Analysis (All Pairs)
# ==========================================================================
run_deseq_by_name_combos <- function(count_matrix, meta_dt, RDS, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, group)])
  rownames(col_data) <- col_data$sample_id
  col_data$group <- as.factor(col_data$group)

  # ✅ FIXED: Uses your Excel "group" column to get clean all-vs-all combinations
  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")

  combos <- combn(levels(col_data$group), 2, simplify = FALSE)

  results_dt <- rbindlist(lapply(combos, function(pair) {
    g2 <- pair[1]
    g1 <- pair[2]
    extract_results(
      dds,
      c("group", g1, g2),
      paste0(g1, "_vs_", g2),
      list(Test_Group = g1, Ref_Group = g2)
    )
  }))
  
  return(list(dt = results_dt, dds = dds))
  return(results_dt)
  message("✓ By group combinations (all combos) : ", nrow(results_dt), " rows generated.")
}

# ==========================================================================
# Execution Core
# ==========================================================================

# 1. Run all analyses and return both the results table AND the dds object from each
# (Make sure your functions return a list(dt = results_dt, dds = dds))
res_ref    <- run_deseq_by_name_ref(count_matrix, meta_dt, REF) 
res_date   <- run_deseq_by_date(count_matrix, meta_dt)
res_combos <- run_deseq_by_name_combos(count_matrix, meta_dt)

# 2. Compile all dds models into a single structured list for the RDS output
master_rds <- list(
  ref   = res_ref$dds,
  date  = res_date$dds,
  combo = res_combos$dds
)
saveRDS(master_rds, RDS)

# 3. Tag and bind all statistical tables together for the Parquet output
dt_ref   <- res_ref$dt[, Contrast_Type := "ref"]
dt_date  <- res_date$dt[, Contrast_Type := "date"]
dt_combo <- res_combos$dt[, Contrast_Type := "combo"]

master_parquet <- rbindlist(list(dt_ref, dt_date, dt_combo), use.names = TRUE, fill = TRUE)
write_parquet(master_parquet, PARQUET)

message("✅ Master RDS and Parquet files compiled successfully.")