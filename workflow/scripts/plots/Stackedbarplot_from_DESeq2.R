################################################################################
# Project : "MicrobExplorer"
# Script: "Stackedbarplot from DESeq2 Results"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(phyloseq)
library(ggplot2)
library(data.table)
library(viridis)
library(arrow)

# Load helpers
source("scripts/plots/utils_Barplot_DESeq2.R")

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Input
DESEQ_FILES   <- snakemake[["input"]][["deseq_files"]]
PHYLOSEQ_OBJ       <- snakemake[["input"]][["phyloseq_obj"]]
METADATA      <- snakemake[["input"]][["metadata"]]

# Output
PDF     <- snakemake[["output"]][["pdf"]]
PARQUET <- snakemake[["output"]][["parquet"]]

# Thresholds from config/params
PADJ_THRESHOLD <- snakemake@params[["padj"]] %||% 0.05
LFC_THRESH     <- snakemake@params[["lfc"]] %||% 0
CONTRASTS_LIST <- snakemake@params[["contrasts"]]  # e.g., ["ref", "date", "combo"]

# ==========================================================================
# 1. Data Loading
# ==========================================================================
message("Loading data...")
message(sprintf("  → DESeq2 contrasts configured: %s", paste(CONTRASTS_LIST, collapse=", ")))

# Load DESeq2 results dynamically from file list
deseq_paths <- list()
for (i in seq_along(DESEQ_FILES)) {
  file_path <- DESEQ_FILES[[i]]
  # Extract contrast type from filename (e.g., "deseq2_ref_contigs.rds" → "ref")
  contrast_name <- gsub(".*deseq2_([^_]+)_.*", "\\1", basename(file_path))
  deseq_paths[[contrast_name]] <- file_path
}

message(sprintf("  → Found %d files: %s", length(deseq_paths), paste(names(deseq_paths), collapse=", ")))

results_list <- load_deseq2_results(deseq_paths)

# Load and prepare phyloseq object
ps <- readRDS(PHYLOSEQ_OBJ)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
df_abundance <- phyloseq_to_dt(ps_rel)

message(sprintf("✓ Loaded phyloseq object: %d OTUs, %d samples", 
                nrow(otu_table(ps)), ncol(otu_table(ps))))

# ==========================================================================
# 2. Generate Plots for Each Contrast
# ==========================================================================
pdf(PDF, width = 14, height = 10)

contrast_names <- names(results_list)

for (contrast_type in contrast_names) {
  
  deseq_result <- results_list[[contrast_type]]
  message(sprintf("\nProcessing %s contrast...", contrast_type))
  
  # Get significant features (up-regulated)
  sig_up <- get_significant_features(deseq_result, PADJ_THRESHOLD, LFC_THRESH, "up")
  
  if (length(sig_up) > 0) {
    message(sprintf("  → %d up-regulated features", length(sig_up)))
    
    # Prepare plot data
    df_plot_up <- prepare_plot_data(df_abundance, sig_up, top_n = 10)
    
    # Create and print plot
    p_up <- create_stackedbarplot(
      df_plot_up,
      title = sprintf("Over-represented Features - %s Contrast", 
                      toupper(contrast_type)),
      subtitle = sprintf("Padj < %g | LogFC > %g", PADJ_THRESHOLD, LFC_THRESH)
    )
    
    print(p_up)
  }
  
  # Get significant features (down-regulated)
  sig_down <- get_significant_features(deseq_result, PADJ_THRESHOLD, LFC_THRESH, "down")
  
  if (length(sig_down) > 0) {
    message(sprintf("  → %d down-regulated features", length(sig_down)))
    
    # Prepare plot data
    df_plot_down <- prepare_plot_data(df_abundance, sig_down, top_n = 10)
    
    # Create and print plot
    p_down <- create_stackedbarplot(
      df_plot_down,
      title = sprintf("Under-represented Features - %s Contrast", 
                      toupper(contrast_type)),
      subtitle = sprintf("Padj < %g | LogFC < -%g", PADJ_THRESHOLD, LFC_THRESH)
    )
    
    print(p_down)
  }
  
  if (length(sig_up) == 0 && length(sig_down) == 0) {
    message(sprintf("  → No significant features found with padj < %g", PADJ_THRESHOLD))
  }
}

dev.off()

message(sprintf("\n✓ PDF generated: %s", PDF))

# ==========================================================================
# 3. Export Results
# ==========================================================================
export_results(df_abundance, PARQUET)

message("\n✓ Processing complete!")