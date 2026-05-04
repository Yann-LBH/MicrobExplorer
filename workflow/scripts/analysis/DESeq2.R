################################################################################
# Project : "MicrobExplorer"
# Script: "DESeq2 apply to Kegg data"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(DESeq2)
library(data.table)
library(arrow)
library(readr)  # uniquement pour read_tsv — fread gère aussi les TSV

# ==========================================================================
# Configuration (variables Snakemake)
# ==========================================================================
data_path     <- snakemake[["input"]][["data_path"]]    # "Statistics/3_Deseq2"
saving_folder <- snakemake[["output"]][["rds_dir"]]     # "RDS_obj/"
parquet_dir   <- snakemake[["parameter"]][["parquet_dir"]] # "Data/Parquet/deseq2/"

dir.create(saving_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(parquet_dir,   recursive = TRUE, showWarnings = FALSE)

setDTthreads(0L)

# ==========================================================================
# 1. Chargement des fichiers — une seule fois pour les 3 analyses
# ==========================================================================
files <- list.files(data_path, pattern = "\\.tsv$", full.names = TRUE)

# Lecture avec data.table (plus rapide que read_tsv en boucle)
# + extraction metadata en une passe
raw_list <- lapply(files, function(f) {
  file_name  <- basename(f)
  digesteur  <- regmatches(file_name, regexpr("TD\\d+", file_name))
  date_raw   <- regmatches(file_name, regexpr("\\d{6}", file_name))
  sample_id  <- paste0(digesteur, "_", date_raw)
  
  dt <- fread(f, showProgress = FALSE)
  setnames(dt, ncol(dt), sample_id)  # renomme la dernière colonne
  
  list(
    dt        = dt,
    sample_id = sample_id,
    digesteur = digesteur,
    date_raw  = date_raw,
    date_obj  = as.Date(date_raw, format = "%y%m%d")
  )
})

# Fusion en une matrice unique — Reduce sur data.table est plus rapide que purrr::reduce
count_data <- Reduce(
  function(a, b) merge(a, b, by = "kegg", all = FALSE),  # inner join
  lapply(raw_list, `[[`, "dt")
)

count_matrix           <- as.matrix(count_data[, -1])
rownames(count_matrix) <- count_data$kegg
count_matrix           <- round(count_matrix)

# Metadata consolidée
meta_dt <- rbindlist(lapply(raw_list, function(x) data.table(
  sample_id = x$sample_id,
  Digesteur = x$digesteur,
  Date_Raw  = x$date_raw,
  Date_Real = x$date_obj
)))
rm(raw_list)

# ==========================================================================
# Fonction utilitaire : extraction des contrastes -> data.table
# ==========================================================================
extract_results <- function(dds, contrast_vec, nom_contraste, extra_cols) {
  res <- results(dds,
                 contrast    = contrast_vec,
                 cooksCutoff = FALSE)
  dt  <- as.data.table(as.data.frame(res), keep.rownames = "KO_Number")
  dt[, Comparison := nom_contraste]
  for (col in names(extra_cols)) dt[, (col) := extra_cols[[col]]]
  dt
}

# ==========================================================================
# 2. Analyse par Digesteur (TD1 comme référence)
# ==========================================================================
run_deseq_by_digesteur_ref <- function(count_matrix, meta_dt, saving_folder, parquet_dir) {
  
  col_data           <- as.data.frame(meta_dt[, .(sample_id, Digesteur)])
  rownames(col_data) <- col_data$sample_id
  col_data$Digesteur <- relevel(as.factor(col_data$Digesteur), ref = "TD1")
  
  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ Digesteur)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  saveRDS(dds, file.path(saving_folder, "diagdds_digesteur_ref.rds"))  # pour PCA/diagnostics
  
  g2     <- "TD1"
  groupes <- levels(col_data$Digesteur)
  
  results_dt <- rbindlist(lapply(
    groupes[groupes != g2],
    function(g1) {
      extract_results(dds,
                      c("Digesteur", g1, g2),
                      paste0(g1, "_vs_", g2),
                      list(Test_Group = g1, Ref_Group = g2))
    }
  ))
  
  write_parquet(results_dt,
                file.path(parquet_dir, "deseq2_digesteur_vs_TD1.parquet"))
  message("✓ Par digesteur (ref TD1) : ", nrow(results_dt), " lignes")
}

# ==========================================================================
# 3. Analyse par date (T vs T-1)
# ==========================================================================
run_deseq_by_date <- function(count_matrix, meta_dt, saving_folder, parquet_dir) {
  
  col_data            <- as.data.frame(meta_dt[, .(sample_id, Date_Raw, Date_Real)])
  rownames(col_data)  <- col_data$sample_id
  col_data$Date_Group <- as.factor(col_data$Date_Raw)
  
  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ Date_Group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  saveRDS(dds, file.path(saving_folder, "diagdds_date.rds"))
  
  # Timeline triée chronologiquement
  timeline <- unique(meta_dt[, .(Date_Real, Date_Raw)])[order(Date_Real)]
  
  if (nrow(timeline) < 2L) {
    warning("Moins de 2 dates : pas de contraste T vs T-1 possible.")
    return(invisible(NULL))
  }
  
  results_dt <- rbindlist(lapply(
    2:nrow(timeline),
    function(i) {
      t_curr <- as.character(timeline$Date_Raw[i])
      t_prev <- as.character(timeline$Date_Raw[i - 1L])
      extract_results(dds,
                      c("Date_Group", t_curr, t_prev),
                      paste0(t_curr, "_vs_", t_prev),
                      list(Date_Test = t_curr, Date_Ref = t_prev))
    }
  ))
  
  write_parquet(results_dt,
                file.path(parquet_dir, "deseq2_date_timeline.parquet"))
  message("✓ Par date (T vs T-1) : ", nrow(results_dt), " lignes")
}

# ==========================================================================
# 4. Analyse par digesteur (toutes combinaisons)
# ==========================================================================
run_deseq_by_digesteur_combos <- function(count_matrix, meta_dt, saving_folder, parquet_dir) {
  
  col_data           <- as.data.frame(meta_dt[, .(sample_id, Digesteur, Date_Raw)])
  rownames(col_data) <- col_data$sample_id
  col_data$Digesteur <- as.factor(col_data$Digesteur)
  
  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~ Digesteur)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  saveRDS(dds, file.path(saving_folder, "diagdds_digesteur_combos.rds"))
  
  combos     <- combn(levels(col_data$Digesteur), 2, simplify = FALSE)
  
  results_dt <- rbindlist(lapply(combos, function(pair) {
    g2 <- pair[1]; g1 <- pair[2]
    extract_results(dds,
                    c("Digesteur", g1, g2),
                    paste0(g1, "_vs_", g2),
                    list(Test_Group = g1, Ref_Group = g2))
  }))
  
  write_parquet(results_dt,
                file.path(parquet_dir, "deseq2_digesteur_combos.parquet"))
  message("✓ Par digesteur (combos) : ", nrow(results_dt), " lignes")
}

# ==========================================================================
# Exécution
# ==========================================================================
run_deseq_by_digesteur_ref   (count_matrix, meta_dt, saving_folder, parquet_dir)
run_deseq_by_date            (count_matrix, meta_dt, saving_folder, parquet_dir)
run_deseq_by_digesteur_combos(count_matrix, meta_dt, saving_folder, parquet_dir)

message("Pipeline DESeq2 terminée.")