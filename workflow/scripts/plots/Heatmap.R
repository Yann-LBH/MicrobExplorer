################################################################################
# Project : "MicrobExplorer"
# Script: "Heatmap"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
suppressPackageStartupMessages({
  # Libraries CRAN
  library(data.table)
  library(rlang)
  library(circlize)
  library(viridis)
  library(vegan)
  library(arrow)
  # Libraries Bioconductor
  library(phyloseq)
  library(ComplexHeatmap)
})

# ==========================================================================
# Configuration Snakemake
# ==========================================================================
source("workflow/scripts/utils/utils_title_resolver.R")

# Inputs
PHYLOSEQ_OBJ <- as.character(snakemake@input[["phyloseq_obj"]])[1]

# Outputs
PDF     <- as.character(snakemake@output[["pdf"]])[1]
PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared plots features
SHARED      <- snakemake@params[["shared"]]
PALETTE     <- as.character(SHARED$palette) %||% "turbo"
PDF_SIZE    <- as.numeric(SHARED$pdf_size) %||% c(14, 12)
TITLE_SIZE  <- as.integer(SHARED$title_size) %||% 12
LEGEND_SIZE <- as.integer(SHARED$legend_size) %||% 10
RANK        <- as.character(snakemake@params[["rank"]])[1]

# Parameters
TOP_N             <- as.integer(snakemake@params[["top_n"]])[1] %||% 50
CLUST_METHOD      <- as.character(snakemake@params[["clust_method"]])[1] %||% "complete"
DISTANCE_METHOD   <- as.character(snakemake@params[["distance_method"]])[1] %||% "bray"

# Wildcards & Variables globales
SOURCE       <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

RESOLVED_TITLE    <- function(sample_id) resolve_text("TITLE_HEATMAP", source = toupper(SOURCE), top_n = TOP_N, sample = sample_id)
RESOLVED_SUBTITLE <- resolve_text("SUBTITLE_HEATMAP", clust_method = CLUST_METHOD, distance_method = DISTANCE_METHOD, rank = RANK)

# ==========================================================================
# 1. Chargement de l'objet Phyloseq
# ==========================================================================
ps <- readRDS(PHYLOSEQ_OBJ)

# ==========================================================================
# 2. Agglomération selon le RANK (Strict : Erreur si RANK invalide)
# ==========================================================================
if (RANK != "") {
  available_ranks <- rank_names(ps)
  
  if (!RANK %in% available_ranks) {
    stop(sprintf(
      "❌ Configuration Error: The requested RANK '%s' was not found in phyloseq object ranks [%s].",
      RANK, paste(available_ranks, collapse = ", ")
    ))
  }
  
  message("--- Aggregating at rank: ", RANK, " ---")
  ps <- tax_glom(ps, taxrank = RANK, NArm = FALSE)
}

# ==========================================================================
# 3. Échelle de couleurs GLOBALE (Calculée après agglomération)
# ==========================================================================
global_mat <- as.matrix(otu_table(ps)@.Data)
if (!taxa_are_rows(ps)) global_mat <- t(global_mat)
global_mat[is.na(global_mat)] <- 0
global_max <- max(log10(global_mat + 1), na.rm = TRUE)
if (is.na(global_max) || global_max <= 0) global_max <- 1

col_fun <- colorRamp2(seq(0, global_max, length.out = 5), viridis(5, option = tolower(PALETTE)))

# ==========================================================================
# 4. Génération Heatmaps (PDF) + Export Tidy Data (Shiny)
# ==========================================================================
samples <- unique(sample_data(ps)$name)
if (anyNA(samples)) {
  warning(sprintf("%d sample(s) with NA 'name' in metadata will be skipped.", sum(is.na(samples))))
  samples <- samples[!is.na(samples)]
}

# Ouverture du PDF avec fermeture automatique sécurisée
pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
on.exit(if (names(dev.cur()) != "null device") dev.off())

shiny_data_list <- lapply(samples, function(s) {
  
  # --- Titles ---
  current_title <- RESOLVED_TITLE(s)
  full_title    <- paste0(current_title, "\n", RESOLVED_SUBTITLE)

  # --- Pruning & Matrix ---
  keep_samples <- sample_names(ps)[sample_data(ps)$name == s]
  ps_sub <- prune_samples(keep_samples, ps)
  
  if (ntaxa(ps_sub) == 0) {
    grid::grid.newpage()
    grid::grid.text(paste("Pas assez de données pour :", s), gp = grid::gpar(fontsize = 14))
    return(NULL)
  }
  
  top_taxa <- names(sort(taxa_sums(ps_sub), decreasing = TRUE))[seq_len(min(TOP_N, ntaxa(ps_sub)))]
  ps_top <- prune_taxa(top_taxa, ps_sub)
  
  mat <- as.matrix(otu_table(ps_top)@.Data)
  if (!taxa_are_rows(ps_top)) mat <- t(mat)
  
  mat[is.na(mat)] <- 0
  mat <- log10(mat + 1)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    grid::grid.newpage()
    grid::grid.text(paste("Pas assez de données pour :", s), gp = grid::gpar(fontsize = 14))
    return(NULL)
  }
  
  # --- Extraction des Labels (Strict selon RANK) ---
  ps_top <- prune_taxa(rownames(mat), ps_top)
  tax_mat_top <- as.matrix(tax_table(ps_top))
  
  target_rank <- if (RANK != "") RANK else colnames(tax_mat_top)[ncol(tax_mat_top)]
  if (!target_rank %in% colnames(tax_mat_top)) {
    stop(sprintf("❌ Configuration Error: Target rank '%s' not present in tax_table columns.", target_rank))
  }
  
  tax_vals <- as.character(tax_mat_top[, target_rank])
  tax_vals[is.na(tax_vals) | tax_vals == "" | tax_vals == "Unassigned"] <- "Unknown"

  # --- 1. Drawing Heatmap in PDF ---
  dist_matrix <- as.matrix(vegdist(mat, method = DISTANCE_METHOD))
  dist_matrix[is.na(dist_matrix)] <- 0
  hc_rows <- hclust(as.dist(dist_matrix), method = CLUST_METHOD)

  draw(
    Heatmap(
      mat,
      column_title = full_title,
      column_title_gp = grid::gpar(fontsize = TITLE_SIZE, fontface = "bold"),
      name = "Abondance\n(log10)",
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = LEGEND_SIZE, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = LEGEND_SIZE),
        legend_height = grid::unit(5, "cm"),
        grid_width = grid::unit(0.8, "cm")
      ),
      cluster_rows = hc_rows,
      cluster_columns = FALSE,
      row_dend_width = grid::unit(25, "mm"),
      row_labels = tax_vals,
      row_names_gp = grid::gpar(fontsize = 9, fontitalic = TRUE),
      column_names_gp = grid::gpar(fontsize = 10),
      col = col_fun
    ),
    padding = grid::unit(c(5, 5, 5, 25), "mm")
  )

  # --- 2. EXTRACTION SECURISEE POUR SHINY ---
  dt_shiny <- as.data.table(as.table(mat))
  setnames(dt_shiny, c("feature_id", "condition_col", "log10_abundance"))

  # Jointure explicite par clé (nom d'entité) au lieu d'un recyclage positionnel
  label_map <- setNames(tax_vals, rownames(mat))

  dt_shiny[, `:=`(
    sample_name   = s,
    feature_label = label_map[as.character(feature_id)],
    raw_abundance = (10^log10_abundance) - 1,
    rank          = target_rank,
    source        = SOURCE
  )]

  return(dt_shiny)
})

# ==========================================================================
# 5. Sauvegarde du fichier Parquet consolidé pour Shiny
# ==========================================================================
heatmap_parquet_dt <- rbindlist(Filter(Negate(is.null), shiny_data_list), fill = TRUE)

if (nrow(heatmap_parquet_dt) == 0) {
  stop("❌ No heatmap could be generated: there was insufficient data.")
}

write_parquet(heatmap_parquet_dt, PARQUET)

message("✓ PDF Multi-pages (échelle globale) : ", PDF)
message("✓ Parquet Shiny                    : ", PARQUET)