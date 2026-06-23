################################################################################
# Project : "MicrobExplorer"
# Script  : "Analysis : DESeq2"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(DESeq2)
library(arrow)
library(readxl)

# Extraction des objets via la syntaxe S4 @
DATA <- snakemake@input[["data"]]
METADATA <- snakemake@input[["metadata"]]
RDS <- snakemake@output[["rds"]]
PARQUET <- snakemake@output[["parquet"]]

# Contrôles et paramètres
CONTRASTS <- snakemake@params[["contrasts"]]
# Si non définie dans Snakemake, prend la première valeur triée de la colonne "name"
REF <- snakemake@params[["ref"]]
# ==========================================================================
# 1. Chargement des Métadonnées et des Fichiers
# ==========================================================================
# Lecture du fichier Excel de métadonnées
meta_dt <- as.data.table(read_xlsx(METADATA))

# Nettoyage et forçage des types de base requis
meta_dt[, sample_id := as.character(sample_id)]
meta_dt[, name      := as.character(name)]
meta_dt[, date      := as.character(date)]
meta_dt[, Date_Real := as.Date(date, format = "%y%m%d")] # Conversion en vraie date

# Gestion dynamique de la référence (REF)
if (is.null(REF) || REF == "" || is.na(REF)) {
  REF <- sort(meta_dt$name)[1]
}

# Liste de tous les fichiers TSV présents dans le dossier d'entrée
files <- list.files(DATA, pattern = "\\.tsv$", full.names = TRUE)

# Lecture croisée dynamique : on charge le fichier s'il matche un sample_id des métadonnées
raw_list <- lapply(files, function(f) {
  file_name <- basename(f)
  
  # On cherche quel sample_id des métadonnées est contenu dans le nom du fichier tsv
  matched_sample <- meta_dt[sapply(sample_id, function(sid) grepl(sid, file_name)), sample_id]
  
  if (length(matched_sample) == 0) return(NULL) # Ignore les fichiers hors métadonnées
  
  dt <- fread(f, showProgress = FALSE)
  # On s'assure que la colonne d'abondance porte le nom du sample_id
  setnames(dt, ncol(dt), matched_sample) 
  
  list(dt = dt, sample_id = matched_sample)
})

# Suppression des éléments NULL
raw_list <- Filter(Negate(is.null), raw_list)

# Fusion des matrices de comptage (Inner join basé sur la colonne "kegg")
count_data <- Reduce(
  function(a, b) merge(a, b, by = "kegg", all = FALSE),
  lapply(raw_list, `[[`, "dt")
)

count_matrix <- as.matrix(count_data[, -1])
rownames(count_matrix) <- count_data$kegg
count_matrix <- round(count_matrix)

# Filtrage final des métadonnées pour ne garder que les échantillons effectivement chargés
meta_dt <- meta_dt[sample_id %in% colnames(count_matrix)]
rm(raw_list)

# ==========================================================================
# Fonction utilitaire : extraction des contrastes
# ==========================================================================
extract_results <- function(dds, contrast_vec, nom_contraste, extra_cols) {
  res <- results(dds, contrast = contrast_vec, cooksCutoff = FALSE)
  dt <- as.data.table(as.data.frame(res), keep.rownames = "KO_Number")
  dt[, Comparison := nom_contraste]
  for (col in names(extra_cols)) dt[, (col) := extra_cols[[col]]]
  dt
}

# ==========================================================================
# 2. Analyse par Groupe d'échantillon (Nom de condition vs Référence)
# ==========================================================================
run_deseq_by_name_ref <- function(count_matrix, meta_dt, REF, RDS, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, name)])
  rownames(col_data) <- col_data$sample_id
  
  # Relevel dynamique basé sur la variable REF
  col_data$name <- relevel(as.factor(col_data$name), ref = REF)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~name)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")
  
  # Sauvegarde du RDS principal pour les analyses en aval
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
  message("✓ Par condition (ref ", REF, ") : ", nrow(results_dt), " lignes")
}

# ==========================================================================
# 3. Analyse par date chronologique (T vs T-1)
# ==========================================================================
run_deseq_by_date <- function(count_matrix, meta_dt, PARQUET) {
  col_data <- as.data.frame(meta_dt[, .(sample_id, date, Date_Real)])
  rownames(col_data) <- col_data$sample_id
  col_data$Date_Group <- as.factor(col_data$date)

  dds <- DESeqDataSetFromMatrix(count_matrix, col_data, design = ~Date_Group)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, test = "Wald", fitType = "parametric")

  # Extraction de la timeline ordonnée
  timeline <- unique(meta_dt[, .(Date_Real, date)])[order(Date_Real)]

  if (nrow(timeline) < 2L) {
    warning("Moins de 2 dates : pas de contraste T vs T-1 possible.")
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

  # Ajout à la suite du fichier Parquet ou écriture séparée selon vos besoins
  # Ici on remplace ou crée un parquet dédié si nécessaire
  write_parquet(results_dt, gsub("\\.parquet$", "_time.parquet", PARQUET))
  message("✓ Par date (T vs T-1) : ", nrow(results_dt), " lignes")
}

# ==========================================================================
# 4. Analyse par combinaisons de noms (Toutes les paires possibles)
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
  message("✓ Par condition (all combos) : ", nrow(results_dt), " lignes")
}

# ==========================================================================
# Exécution du Pipeline Réorganisé
# ==========================================================================
run_deseq_by_name_ref(count_matrix, meta_dt, REF, RDS, PARQUET)
run_deseq_by_date(count_matrix, meta_dt, PARQUET)
run_deseq_by_name_combos(count_matrix, meta_dt, PARQUET)

message("✓ Pipeline DESeq2 terminée avec succès.")