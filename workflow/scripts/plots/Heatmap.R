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
try(snakemake, silent = TRUE)
if (!exists("snakemake")) {
  snakemake <- list(
    input  = list(
      rds  = "Data/RDS/phyloseq.rds",
      tsv  = "Data/heatmap_data.tsv"
    ),
    output = list(
      pdf     = "Graphique/Heatmap/Heatmaps.pdf",
      parquet = "Data/Parquet/heatmap_data.parquet"
    ),
    params = list(
      top_n        = 50L,
      color_opt    = "turbo",
      clust_method = "ward.D2"
    )
  )
}

setDTthreads(0L)

top_n        <- as.integer(snakemake$params$top_n)
color_opt    <- snakemake$params$color_opt
clust_method <- snakemake$params$clust_method

dir.create(dirname(snakemake$output$pdf),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake$output$parquet), recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 1. Chargement
# ==========================================================================
ps      <- readRDS(snakemake$input$rds)
dt      <- fread(snakemake$input$tsv, sep = "\t")

# ==========================================================================
# 2. Export Parquet pour Shiny
# ==========================================================================
write_parquet(dt, snakemake$output$parquet)

# ==========================================================================
# 3. Heatmaps -> PDF compilé
# ==========================================================================
reacteurs <- unique(sample_data(ps)$Digesteur)

pdf(snakemake$output$pdf, width = 14, height = 12)

lapply(reacteurs, function(r) {
  
  ps_sub   <- subset_samples(ps, Digesteur == r)
  top_taxa <- names(sort(taxa_sums(ps_sub), decreasing = TRUE))[seq_len(min(top_n, ntaxa(ps_sub)))]
  ps_top   <- prune_taxa(top_taxa, ps_sub)
  
  mat <- log10(otu_table(ps_top)@.Data + 1)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  if (nrow(mat) < 2L) {
    grid::grid.newpage()
    grid::grid.text(paste("Pas assez de données pour", r))
    return(invisible(NULL))
  }
  
  ps_top    <- prune_taxa(rownames(mat), ps_top)
  dist_bray <- vegdist(mat, method = "bray")
  dist_bray[is.na(dist_bray)] <- 0
  
  col_fun <- colorRamp2(
    seq(0, max(mat, na.rm = TRUE), length.out = 5),
    viridis(5, option = color_opt)
  )
  
  draw(
    Heatmap(
      mat,
      column_title         = paste("Réacteur :", r),
      column_title_gp      = grid::gpar(fontface = "bold"),
      name                 = "Abondance (log10)",
      heatmap_legend_param = list(
        title_gp      = grid::gpar(fontsize = 12, fontface = "bold"),
        labels_gp     = grid::gpar(fontsize = 10),
        legend_height = grid::unit(6, "cm"),
        grid_width    = grid::unit(1, "cm")
      ),
      clustering_distance_rows = dist_bray,
      clustering_method_rows   = clust_method,
      cluster_columns          = FALSE,
      row_dend_width           = grid::unit(20, "mm"),
      row_labels               = as.character(tax_table(ps_top)[, "species"]),
      row_names_gp             = grid::gpar(fontsize = 10),
      col                      = col_fun
    ),
    padding = grid::unit(c(2, 2, 2, 20), "mm")
  )
})

dev.off()

message("✓ PDF     : ", snakemake$output$pdf)
message("✓ Parquet : ", snakemake$output$parquet)