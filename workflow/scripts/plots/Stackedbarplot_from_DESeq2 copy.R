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
# Input files
deseq_ref_path  <- snakemake[["input"]][["deseq_ref"]]
deseq_date_path <- snakemake[["input"]][["deseq_date"]]
deseq_combo_path <- snakemake[["input"]][["deseq_combo"]]
ps_path         <- snakemake[["input"]][["phyloseq"]]

# Output files
pdf_out     <- snakemake[["output"]][["pdf"]]
parquet_out <- snakemake[["output"]][["parquet"]]

# Thresholds from config/params
padj_threshold <- snakemake@params[["padj"]] %||% 0.05
lfc_thresh     <- snakemake@params[["lfc"]] %||% 0

# ==========================================================================
# 1. Data Loading
# ==========================================================================
message("Loading data...")

# Load DESeq2 results
deseq_paths <- list(
  deseq_ref   = deseq_ref_path,
  deseq_date  = deseq_date_path,
  deseq_combo = deseq_combo_path
)

results_list <- load_deseq2_results(deseq_paths)

# Load and prepare phyloseq object
ps <- readRDS(ps_path)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
df_abundance <- phyloseq_to_dt(ps_rel)

message(sprintf("✓ Loaded phyloseq object: %d OTUs, %d samples", 
                nrow(otu_table(ps)), ncol(otu_table(ps))))

# ==========================================================================
# 2. Generate Plots for Each Contrast
# ==========================================================================
pdf(pdf_out, width = 14, height = 10)

contrast_names <- names(results_list)

for (contrast_type in contrast_names) {
  
  deseq_result <- results_list[[contrast_type]]
  message(sprintf("\nProcessing %s contrast...", contrast_type))
  
  # Get significant features (up-regulated)
  sig_up <- get_significant_features(deseq_result, padj_threshold, lfc_thresh, "up")
  
  if (length(sig_up) > 0) {
    message(sprintf("  → %d up-regulated features", length(sig_up)))
    
    # Prepare plot data
    df_plot_up <- prepare_plot_data(df_abundance, sig_up, top_n = 10)
    
    # Create and print plot
    p_up <- create_stackedbarplot(
      df_plot_up,
      title = sprintf("Over-represented Features - %s Contrast", 
                      toupper(contrast_type)),
      subtitle = sprintf("Padj < %g | LogFC > %g", padj_threshold, lfc_thresh)
    )
    
    print(p_up)
  }
  
  # Get significant features (down-regulated)
  sig_down <- get_significant_features(deseq_result, padj_threshold, lfc_thresh, "down")
  
  if (length(sig_down) > 0) {
    message(sprintf("  → %d down-regulated features", length(sig_down)))
    
    # Prepare plot data
    df_plot_down <- prepare_plot_data(df_abundance, sig_down, top_n = 10)
    
    # Create and print plot
    p_down <- create_stackedbarplot(
      df_plot_down,
      title = sprintf("Under-represented Features - %s Contrast", 
                      toupper(contrast_type)),
      subtitle = sprintf("Padj < %g | LogFC < -%g", padj_threshold, lfc_thresh)
    )
    
    print(p_down)
  }
  
  if (length(sig_up) == 0 && length(sig_down) == 0) {
    message(sprintf("  → No significant features found with padj < %g", padj_threshold))
  }
}

dev.off()

message(sprintf("\n✓ PDF generated: %s", pdf_out))

# ==========================================================================
# 3. Export Results
# ==========================================================================
export_results(df_abundance, parquet_out)

message("\n✓ Processing complete!")