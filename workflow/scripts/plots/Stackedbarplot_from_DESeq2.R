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
library(rlang)

# Libraries Bioconductor
library(phyloseq)

# Load helper
helper_path <- dirname(snakemake@scriptdir)
source(file.path(helper_path, "utils", "utils_Barplot_DESeq2.R"))

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Input
DESEQ_FILES <- as.character(c(snakemake@input[["deseq_files"]])[1])
PHYLOSEQ_OBJ <- as.character(c(snakemake@input[["phyloseq_obj"]])[1])
METADATA <- as.character(c(snakemake@input[["metadata"]])[1])

# Output
PDF <- as.character(c(snakemake@output[["pdf"]])[1])
PARQUET <- as.character(c(snakemake@output[["parquet"]])[1])

# Thresholds from config/params
TITLE <- as.character(c(snakemake@params[["title"]])[1]) %||% "Top N Features Differentially Abundant"
SUBTITLE <- as.character(c(snakemake@params[["subtitle"]])[1]) %||% "DESeq2 Results"
PADJ_THRESHOLD <- as.numeric(c(snakemake@params[["padj"]])[1]) %||% 0.05
LFC_THRESH <- as.numeric(c(snakemake@params[["lfc"]])[1]) %||% 0
CONTRASTS_LIST <- unlist(as.character(c(snakemake@params[["contrasts"]])[1])) # e.g., ["ref", "date", "combo"]
TOP_N <- as.integer(c(snakemake@params[["top_n"]])[1]) %||% 10

# ==========================================================================
# 1. Data Loading
# ==========================================================================
message("Loading data...")

master_table <- as.data.table(read_parquet(DESEQ_FILES))

results_list <- split(master_table, master_table$Contrast_Type)

message(sprintf("  → Successfully split master parquet into %d contrast views: %s", 
                length(results_list), paste(names(results_list), collapse = ", ")))

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

  # Extract significant up-regulated features
  sig_up <- get_significant_features(deseq_result, PADJ_THRESHOLD, LFC_THRESH, "up")

  if (length(sig_up) > 0) {
    message(sprintf("  → %d up-regulated features", length(sig_up)))

    df_plot_up <- prepare_plot_data(df_abundance, sig_up, TOP_N = TOP_N)

    resolved_title <- sprintf(TITLE, toupper(contrast_type))
    resolved_subtitle <- sprintf(SUBTITLE, PADJ_THRESHOLD, LFC_THRESH)

    p_up <- create_stackedbarplot(
      df_plot_up,
      title = resolved_title,
      subtitle = resolved_subtitle
    )
    print(p_up)
  }

  # Extract significant down-regulated features
  sig_down <- get_significant_features(deseq_result, PADJ_THRESHOLD, LFC_THRESH, "down")

  if (length(sig_down) > 0) {
    message(sprintf("  → %d down-regulated features", length(sig_down)))

    df_plot_down <- prepare_plot_data(df_abundance, sig_down, TOP_N = TOP_N)

    resolved_title_down <- sprintf(TITLE, toupper(contrast_type))
    resolved_subtitle_down <- sprintf(SUBTITLE, PADJ_THRESHOLD, -LFC_THRESH)

    p_down <- create_stackedbarplot(
      df_plot_down,
      title    = resolved_title_down,
      subtitle = resolved_subtitle_down
    )
    print(p_down)
  }

  if (length(sig_up) == 0 && length(sig_down) == 0) {
    message(sprintf("  → No significant features found with padj < %g", PADJ_THRESHOLD))
  }
}

dev.off()
message(sprintf("\n✓ PDF : Stackedbarplot from DESeq2 Successfully generated: %s", PDF))

# ==========================================================================
# 3. Export Results
# ==========================================================================
export_results(df_abundance, PARQUET)
message("✓ PARQUET : Stackedbarplot from DESeq2 Successfully generated.")