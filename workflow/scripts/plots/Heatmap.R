################################################################################
# Project : "MicrobExplorer"
# Script: "Heatmap"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(rlang)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(vegan)
library(arrow)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================

# Inputs
DATA <- as.character(snakemake@input[["data"]])
PHYLOSEQ_OBJ <- as.character(c(snakemake@input[["phyloseq_obj"]])[1])
METADATA <- as.character(c(snakemake@input[["metadata"]])[1])

# Outputs
PDF <- as.character(c(snakemake@output[["pdf"]])[1])
PARQUET <- as.character(c(snakemake@output[["parquet"]])[1])

# Parameters
TOP_N <- as.integer(c(snakemake@params[["top_n"]])[1])
TAXON_RANK <- as.character(snakemake@params[["taxon_rank"]][1])
COLOR_OPT <- as.character(c(snakemake@params[["color_opt"]])[1]) %||% "turbo"
CLUST_METHOD <- as.character(c(snakemake@params[["clust_method"]])[1]) %||% "complete"
DISTANCE_METHOD <- as.character(c(snakemake@params[["distance_method"]])[1]) %||% "bray"


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
# 2. Heatmaps Generation -> Multi-page PDF
# ==========================================================================

# FIXED: Standardizing column reference to 'name' as seen in metadata Excel (e.g. TD1, TD2)
if (!"name" %in% colnames(sample_data(ps))) {
  # Fallback to Digesteur if 'name' was not mapped during phyloseq object creation
  if ("Digesteur" %in% colnames(sample_data(ps))) {
    sample_data(ps)$name <- sample_data(ps)$Digesteur
  } else {
    stop("Could not find reactor column ('name' or 'Digesteur') in phyloseq sample_data.")
  }
}

reacteurs <- unique(sample_data(ps)$name)

# Ensure output directory exists
dir.create(dirname(PDF), recursive = TRUE, showWarnings = FALSE)

# Open PDF device safely
pdf(PDF, width = 14, height = 12)

lapply(reacteurs, function(r) {
  # Subset phyloseq for the current reactor
  keep_samples <- sample_names(ps)[sample_data(ps)$name == r]
  ps_sub <- prune_samples(keep_samples, ps)
  
  if (ntaxa(ps_sub) == 0) return(invisible(NULL))
  
  # Select Top N taxa for this specific reactor
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
    grid::grid.text(paste("Pas assez de variations/données pour le réacteur :", r), gp = grid::gpar(fontsize = 14))
    return(invisible(NULL))
  }
  
  # Re-align tax table with remaining matrix rows
  ps_top <- prune_taxa(rownames(mat), ps_top)
  
  # Extract species taxonomy labels cleanly, replace NA with Unknown
  tax_labels <- as.character(tax_table(ps_top)[, TAXON_RANK])
  
  tax_labels <- paste0(rownames(mat), " (", tax_labels, ")")

  # FIXED: Compute a clean, robust Bray-Curtis distance matrix for rows
  dist_matrix <- as.matrix(vegdist(mat, method = DISTANCE_METHOD))
  dist_matrix[is.na(dist_matrix)] <- 0
  hc_rows <- hclust(as.dist(dist_matrix), method = CLUST_METHOD)
  
  # Dynamic Color Mapping
  max_val <- max(mat, na.rm = TRUE)
  if (is.na(max_val) || max_val <= 0) {
    max_val <- 1 # Fallback value to avoid seq(0, 0) crash
  }

  # Create a clean, evaluated numeric vector for color mapping
  col_fun <- colorRamp2(
    seq(0, max_val, length.out = 5),
    viridis(5, option = COLOR_OPT)
  )
  
  # Draw Heatmap onto PDF page
  draw(
    Heatmap(
      mat,
      column_title = paste("Réacteur :", r),
      column_title_gp = grid::gpar(fontsize = 16, fontface = "bold"),
      name = "Abondance\n(log10)",
      heatmap_legend_param = list(
        title_gp      = grid::gpar(fontsize = 11, fontface = "bold"),
        labels_gp     = grid::gpar(fontsize = 9),
        legend_height = grid::unit(5, "cm"),
        grid_width    = grid::unit(0.8, "cm")
      ),
      cluster_rows = hc_rows,
      cluster_columns = FALSE,
      row_dend_width = grid::unit(25, "mm"),
      row_labels = tax_labels,
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