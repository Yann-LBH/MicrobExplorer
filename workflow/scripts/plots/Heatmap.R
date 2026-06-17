################################################################################
# Project : "MicrobExplorer"
# Script: "Heatmap"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
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
DATA <- as.character(snakemake@input$data)
PHYLOSEQ_OBJ <- as.character(snakemake@input$phyloseq_obj)
METADATA <- as.character(snakemake@input$metadata)

# Outputs
PDF <- as.character(snakemake@output$pdf)
PARQUET <- as.character(snakemake@output$parquet)

# Parameters
TOP_N <- as.integer(snakemake@params$top_n)
COLOR_OPT <- as.character(snakemake@params$color_opt)
CLUST_METHOD <- as.character(snakemake@params$clust_method)


# ==========================================================================
# 1. Chargement
# ==========================================================================
ps <- readRDS(PHYLOSEQ_OBJ)
dt <- fread(DATA, sep = "\t")

# ==========================================================================
# 2. Export Parquet pour Shiny
# ==========================================================================
write_parquet(dt, PARQUET)

# ==========================================================================
# 3. Heatmaps -> PDF compilé
# ==========================================================================
reacteurs <- unique(sample_data(ps)$Digesteur)

pdf(PDF, width = 14, height = 12)

lapply(reacteurs, function(r) {
  ps_sub <- subset_samples(ps, Digesteur == r)
  top_taxa <- names(sort(taxa_sums(ps_sub), decreasing = TRUE))[seq_len(min(TOP_N, ntaxa(ps_sub)))]
  ps_top <- prune_taxa(top_taxa, ps_sub)

  mat <- log10(otu_table(ps_top)@.Data + 1)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  if (nrow(mat) < 2L) {
    grid::grid.newpage()
    grid::grid.text(paste("Pas assez de données pour", r))
    return(invisible(NULL))
  }

  ps_top <- prune_taxa(rownames(mat), ps_top)
  dist_bray <- vegdist(mat, method = "bray")
  dist_bray[is.na(dist_bray)] <- 0

  col_fun <- colorRamp2(
    seq(0, max(mat, na.rm = TRUE), length.out = 5),
    viridis(5, option = COLOR_OPT)
  )

  draw(
    Heatmap(
      mat,
      column_title = paste("Réacteur :", r),
      column_title_gp = grid::gpar(fontface = "bold"),
      name = "Abondance (log10)",
      heatmap_legend_param = list(
        title_gp      = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp     = grid::gpar(fontsize = 10),
        legend_height = grid::unit(6, "cm"),
        grid_width    = grid::unit(1, "cm")
      ),
      clustering_distance_rows = dist_bray,
      clustering_method_rows = CLUST_METHOD,
      cluster_columns = FALSE,
      row_dend_width = grid::unit(20, "mm"),
      row_labels = as.character(tax_table(ps_top)[, "species"]),
      row_names_gp = grid::gpar(fontsize = 10),
      col = col_fun
    ),
    padding = grid::unit(c(2, 2, 2, 20), "mm")
  )
})

dev.off()

message("✓ PDF     : ", PDF)
message("✓ Parquet : ", PARQUET)
