################################################################################
# Project : "MicrobExplorer"
# Script  : "Analysis : DESeq2"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(readxl)
library(arrow)
library(DESeq2)

DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])
RDS      <- as.character(snakemake@output[["rds"]])
PARQUET  <- as.character(snakemake@output[["parquet"]])

# Controls and parameters
CONTRASTS <- tolower(as.character(snakemake@params[["contrasts"]]))[1]
REF       <- tolower(as.character(snakemake@params[["ref"]]))[1]

# ==========================================================================
# 1. Loading Metadata and Files
# ==========================================================================
meta_dt <- as.data.table(read_xlsx(METADATA))

# ✅ FIXED: Corrected parsing format to match the "dd/mm/yyyy" layout from Excel
meta_dt[, Date_Real := as.Date(date, format = "%d/%m/%Y")]

# Handle dynamic reference defaults
if (is.null(REF) || REF == "" || is.na(REF)) {
  REF <- sort(meta_dt$name)[1]
}

# Scan input TSV files
files <- list.files(DATA, pattern = "\\.tsv$", full.names = TRUE)

# Cross-load files matching existing sample ids
raw_list <- lapply(files, function(f) {
  file_name <- basename(f)
  
  matched_sample <- meta_dt[sapply(sample_id, function(sid) grepl(sid, file_name)), sample_id]
  if (length(matched_sample) == 0) return(NULL) 
  
  dt <- fread(f, showProgress = FALSE)
  
  # ✅ FIXED: Force lowercase column headers to match Python step 6 updates
  setnames(dt, tolower(names(dt)))
  setnames(dt, ncol(dt), matched_sample) 
  
  list(dt = dt, sample_id = matched_sample)
})

raw_list <- Filter(Negate(is.null), raw_list)

# ✅ FIXED: Dynamically find the ID column (contig_id, read_id, or ko) instead of hardcoding "kegg"
sample_headers <- unique(unlist(lapply(raw_list, `[[`, "dt")))
id_col_candidate <- intersect(c("contig_id", "read_id", "ko", "kegg"), names(raw_list[[1]]$dt))[1]

# Merge count matrices using the discovered key identifier
count_data <- Reduce(
  function(a, b) merge(a, b, by = id_col_candidate, all = FALSE),
  lapply(raw_list, `[[`, "dt")
)

count_matrix <- as.matrix(count_data[, !id_col_candidate, with = FALSE])
rownames(count_matrix) <- count_data[[id_col_candidate]]
count_matrix <- round(count_matrix)

# Synchronize metadata structure with loaded matrix values
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
# 2. Condition Analysis (Condition Name vs Reference)
# ==========================================================================
run_deseq_by_name_ref <- function(count_matrix, meta_dt, REF, RDS, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, name)])
  rownames(col_data) <- col_data$sample_id
  
  col_data$name <- relevel(as.factor(col_data$name), ref = REF)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~name)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  
  saveRDS(dds, RDS) 

  groupes <- levels(col_data$name)

  results_dt <- rbindlist(lapply(
    groupes[groupes != REF],
    function(g1) {
      extract_results(
        dds,
        c("name", g1, REF),
        paste0(g1, "_vs_", REF),
        list(Test_Group = g1, Ref_Group = REF)
      )
    }
  ))

  write_parquet(results_dt, PARQUET)
  message("✓ By condition (ref ", REF, ") : ", nrow(results_dt), " rows generated.")
}

# ==========================================================================
# 3. Chronological Analysis (T vs T-1)
# ==========================================================================
run_deseq_by_date <- function(count_matrix, meta_dt, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, date, Date_Real)])
  rownames(col_data) <- col_data$sample_id
  col_data$Date_Group <- as.factor(col_data$date)

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

  write_parquet(results_dt, gsub("\\.parquet$", "_time.parquet", PARQUET))
  message("✓ By timeline (T vs T-1) : ", nrow(results_dt), " rows generated.")
}

# ==========================================================================
# 4. Pairwise Combination Analysis (All Pairs)
# ==========================================================================
run_deseq_by_name_combos <- function(count_matrix, meta_dt, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, name)])
  rownames(col_data) <- col_data$sample_id
  col_data$name <- as.factor(col_data$name)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~name)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")

  combos <- combn(levels(col_data$name), 2, simplify = FALSE)

  results_dt <- rbindlist(lapply(combos, function(pair) {
    g2 <- pair[1]
    g1 <- pair[2]
    extract_results(
      dds,
      c("name", g1, g2),
      paste0(g1, "_vs_", g2),
      list(Test_Group = g1, Ref_Group = g2)
    )
  }))

  write_parquet(results_dt, gsub("\\.parquet$", "_combos.parquet", PARQUET))
  message("✓ By group combinations (all combos) : ", nrow(results_dt), " rows generated.")
}

# ==========================================================================
# Execution Core
# ==========================================================================
run_deseq_by_name_ref(count_matrix, meta_dt, REF, RDS, PARQUET)
run_deseq_by_date(count_matrix, meta_dt, PARQUET)
run_deseq_by_name_combos(count_matrix, meta_dt, PARQUET)

message("✓ DESeq2 pipeline workflow finished successfully.")