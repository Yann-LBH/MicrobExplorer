################################################################################
# Project : "MicrobExplorer"
# Script: "Heatmap"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
library(data.table)
library(rlang)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(vegan)
library(arrow)
library(glue)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================

# Inputs
DATA         <- as.character(snakemake@input[["data"]])
PHYLOSEQ_OBJ <- as.character(snakemake@input[["phyloseq_obj"]])[1]
METADATA     <- as.character(snakemake@input[["metadata"]])[1]

# Outputs
PDF     <- as.character(snakemake@output[["pdf"]])[1]
PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared plots features
SHARED      <- snakemake@params[["shared"]]
PALETTE     <- as.character(SHARED$palette) %||% "turbo"
PDF_SIZE    <- as.numeric(SHARED$pdf_size) %||% c(14, 12)
TITLE_SIZE  <- as.integer(SHARED$title_size) %||% 12
SUBTITLE_SIZE <- as.integer(SHARED$subtitle_size) %||% 10
LEGEND_SIZE <- as.integer(SHARED$legend_size) %||% 10
RANK        <- as.character(snakemake@params[["rank"]])[1]

# Parameters
TITLE_TEMPLATE    <- as.character(snakemake@params[["title"]])[1] %||% "{source} | Abundance of the top {top_n} taxons in sample {sample}"
SUBTITLE_TEMPLATE <- as.character(snakemake@params[["subtitle"]])[1] %||% "Hierarchical clustering: {clust_method} | Distance metric: {distance_method} | {rank}"
TOP_N             <- as.integer(snakemake@params[["top_n"]])[1] %||% 50
CLUST_METHOD      <- as.character(snakemake@params[["clust_method"]])[1] %||% "complete"
DISTANCE_METHOD   <- as.character(snakemake@params[["distance_method"]])[1] %||% "bray"

# Wildcards
SOURCE  <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# Harmonisation de l'affichage du type d'entités
FEATURE_TYPE <- if (SOURCE == "kegg") "pathways" else SOURCE

# ==========================================================================
# 1. Chargement
# ==========================================================================
ps <- readRDS(PHYLOSEQ_OBJ)
tsv_list <- lapply(DATA, function(f) {
  if (file.exists(f)) {
    dt_file <- fread(f, sep = "\t")
    
    # Track the file origin to preserve longitudinal and reactor context
    dt_file[, Sample_Source := basename(f)]
    return(dt_file)
  }
  return(NULL)
})

# Combine everything (keeps all rows, columns, and timepoints dynamically)
annotated_data <- rbindlist(tsv_list, use.names = TRUE, fill = TRUE)

# Export the complete tracking dataset to Parquet
write_parquet(annotated_data, PARQUET)

# ==========================================================================
# 2. Aggregating at designated RANK (Prevents multi-contig saturation)
# ==========================================================================
if (RANK != "" && RANK %in% rank_names(ps)) {
  message("--- Aggregating phyloseq object at rank: ", RANK, " ---")
  # NA taxa elements are preserved or cleaned down by tax_glom
  ps <- tax_glom(ps, taxrank = RANK, NArm = FALSE)
}

# ==========================================================================
# 3. Heatmaps Generation -> Multi-page PDF
# ==========================================================================

samples <- unique(sample_data(ps)$name)

# Open PDF device dynamically sized
pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])

lapply(samples, function(s) {
  # Subset phyloseq for the current sample
  keep_samples <- sample_names(ps)[sample_data(ps)$name == s]
  ps_sub <- prune_samples(keep_samples, ps)
  
  if (ntaxa(ps_sub) == 0) return(invisible(NULL))
  
  # Select Top N features/taxa for this specific reactor *after* aggregation
  top_taxa <- names(sort(taxa_sums(ps_sub), decreasing = TRUE))[seq_len(min(TOP_N, ntaxa(ps_sub)))]
  ps_top <- prune_taxa(top_taxa, ps_sub)
  
  # Extract matrix from OTU table
  mat <- as.matrix(otu_table(ps_top)@.Data)
  
  # CRITICAL FIXED: Ensure Taxa are always ROWS, and Samples are always COLUMNS
  if (!taxa_are_rows(ps_top)) {
    mat <- t(mat)
  }
  
  # Log transformation (log10 + 1 pseudo-count)
  mat <- log10(mat + 1)
  
  # Remove rows (taxa) with zero variance or zero total counts across these subsetted samples
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  # Security check: ComplexHeatmap requires at least 2 rows and 2 columns to perform clustering
  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    grid::grid.newpage()
    grid::grid.text(paste("Pas assez de variations/données pour le sample :", s), gp = grid::gpar(fontsize = 14))
    return(invisible(NULL))
  }
  
  # Re-align tax table with remaining matrix rows
  ps_top <- prune_taxa(rownames(mat), ps_top)
  
  # Extract taxonomy labels cleanly, replace NA with Unknown
  tax_vals <- as.character(tax_table(ps_top)[, RANK])
  tax_vals[is.na(tax_vals) | tax_vals == ""] <- "Unknown"

  # FIXED: Compute a clean, robust distance matrix for rows
  dist_matrix <- as.matrix(vegdist(mat, method = DISTANCE_METHOD))
  dist_matrix[is.na(dist_matrix)] <- 0
  hc_rows <- hclust(as.dist(dist_matrix), method = CLUST_METHOD)
  
  # Dynamic Color Mapping
  max_val <- max(mat, na.rm = TRUE)
  if (is.na(max_val) || max_val <= 0) {
    max_val <- 1 
  }

  col_fun <- colorRamp2(
    seq(0, max_val, length.out = 5),
    viridis(5, option = tolower(PALETTE))
  )
  
  # Resolving layout templates with grid templates
  resolved_title    <- glue(TITLE_TEMPLATE, source = toupper(SOURCE), top_n = TOP_N, sample = samples)
  resolved_subtitle <- glue(SUBTITLE_TEMPLATE, clust_method = CLUST_METHOD, distance_method = DISTANCE_METHOD, rank = RANK)
  full_title <- paste0(resolved_title, "\n", resolved_subtitle)

  # Draw Heatmap onto PDF page
  draw(
    Heatmap(
      mat,
      column_title = full_title,
      column_title_gp = grid::gpar(fontsize = TITLE_SIZE, fontface = "bold"),
      #subtitle_gp = grid::gpar(fontsize = SUBTITLE_SIZE, fontitalic = TRUE),
      name = "Abondance\n(log10)",
      heatmap_legend_param = list(
        title_gp      = grid::gpar(fontsize = LEGEND_SIZE, fontface = "bold"),
        labels_gp     = grid::gpar(fontsize = LEGEND_SIZE),
        legend_height = grid::unit(5, "cm"),
        grid_width    = grid::unit(0.8, "cm")
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
})

dev.off()

message("✓ PDF     : ", PDF)
message("✓ Parquet : ", PARQUET)