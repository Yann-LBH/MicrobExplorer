################################################################################
# Project : "MicrobExplorer"
# Script: "Stackedbarplot from DESeq2 Results"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
library(data.table)
library(ggplot2)
library(viridis)
library(arrow)
library(purrr)

# Libraries Bioconductor
library(phyloseq)

# Load helper
helper_path <- dirname(snakemake@scriptdir)
source(file.path(helper_path, "utils", "utils_Barplot_DESeq2.R"))

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Input
DESEQ_FILES <- as.character(snakemake@input$deseq_files)
PHYLOSEQ_OBJ <- as.character(snakemake@input$phyloseq_obj)
METADATA <- as.character(snakemake@input$metadata)

# Output
PDF <- as.character(snakemake@output$pdf)
PARQUET <- as.character(snakemake@output$parquet)

# Thresholds from config/params
TITLE <- as.character(snakemake@params$title) %||% "Top N Features Differentially Abundant"
SUBTITLE <- as.character(snakemake@params$subtitle) %||% "DESeq2 Results"
PADJ_THRESHOLD <- as.numeric(snakemake@params$padj) %||% 0.05
LFC_THRESH <- as.numeric(snakemake@params$lfc) %||% 0
CONTRASTS_LIST <- unlist(as.character(snakemake@params$contrasts)) # e.g., ["ref", "date", "combo"]
TOP_N <- as.integer(snakemake@params$top_n) %||% 10

# ==========================================================================
# 1. Data Loading
# ==========================================================================
message("Loading data...")
message(sprintf("  → DESeq2 contrasts configured: %s", paste(CONTRASTS_LIST, collapse = ", ")))

# Load DESeq2 results dynamically from file list
deseq_paths <- list()
files_vector <- unlist(DESEQ_FILES)

deseq_paths <- lapply(files_vector, function(file_path) {
  return(file_path)
})

# Extract contrast names systematically using regex matching on filenames
contrast_keys <- gsub(".*deseq2_([^_]+)_.*", "\\1", basename(files_vector))
names(deseq_paths) <- contrast_keys

message(sprintf("  → Found %d files: %s", length(deseq_paths), paste(names(deseq_paths), collapse = ", ")))

results_list <- load_deseq2_results(deseq_paths)

# Load and prepare phyloseq object
ps <- readRDS(PHYLOSEQ_OBJ)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
df_abundance <- phyloseq_to_dt(ps_rel)

message(sprintf(
  "✓ Loaded phyloseq object: %d OTUs, %d samples",
  nrow(otu_table(ps)), ncol(otu_table(ps))
))

# ==========================================================================
# 2. Generate Plots for Each Contrast
# ==========================================================================
pdf(PDF, width = 14, height = 10)

contrast_names <- names(results_list)

for (contrast_type in contrast_names) {
  deseq_result <- results_list[[contrast_type]]
  message(sprintf("\nProcessing %s contrast...", contrast_type))

  # linked to utils_Barplot_DESeq2.R function to prepare data for plotting
  sig_up <- get_significant_features(deseq_result, PADJ_THRESHOLD, LFC_THRESH, "up")

  if (length(sig_up) > 0) {
    message(sprintf("  → %d up-regulated features", length(sig_up)))

    # linked to utils_Barplot_DESeq2.R function to prepare data for plotting
    df_plot_up <- prepare_plot_data(df_abundance, sig_up, TOP_N = TOP_N)

    # DYNAMIC PARSING: Combine contrast tags and configuration keys into labels
    resolved_title <- sprintf(CUSTOM_TITLE, toupper(contrast_type))
    resolved_subtitle <- sprintf(CUSTOM_SUBTITLE, PADJ_THRESHOLD, LFC_THRESH)

    # Create and print plot
    p_up <- create_stackedbarplot(
      df_plot_up,
      title = resolved_title,
      subtitle = resolved_subtitle,
      TOP_N = TOP_N
    )

    print(p_up)
  }

  # Get significant features (down-regulated)
  sig_down <- get_significant_features(deseq_result, PADJ_THRESHOLD, LFC_THRESH, "down")

  if (length(sig_down) > 0) {
    message(sprintf("  → %d down-regulated features", length(sig_down)))

    # linked to utils_Barplot_DESeq2.R function to prepare data for plotting
    df_plot_down <- prepare_plot_data(df_abundance, sig_down, TOP_N = TOP_N)

    # DYNAMIC PARSING: Combine contrast tags and configuration keys into labels
    resolved_title_down <- sprintf(CUSTOM_TITLE, toupper(contrast_type))
    resolved_subtitle_down <- sprintf(CUSTOM_SUBTITLE, PADJ_THRESHOLD, -LFC_THRESH)

    p_down <- create_stackedbarplot(
      df_plot_down,
      title        = resolved_title_down,
      subtitle     = resolved_subtitle_down,
      TOP_N        = TOP_N
    )
    print(p_down)
  }

  if (length(sig_up) == 0 && length(sig_down) == 0) {
    message(sprintf("  → No significant features found with padj < %g", PADJ_THRESHOLD))
  }
}

dev.off()
message(sprintf("\n✓ PDF : Stackedbarplot from DESeq2 Successfully generated.", PDF))

# ==========================================================================
# 3. Export Results
# ==========================================================================
export_results(df_abundance, PARQUET)
message("\n✓ PARQUET : Stackedbarplot from DESeq2 Successfully generated.")
