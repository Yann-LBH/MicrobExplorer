################################################################################
# Project : "MicrobExplorer"
# Script: "Volcano plot from Kegg data DESeq2"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
library(data.table)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(arrow)
library(rlang)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Inputs
DESEQ_FILES <- as.character(snakemake@input[["deseq_files"]])

# Outputs
PARQUET <- as.character(snakemake@output[["parquet"]])[1] # Single parquet file path
PDF <- as.character(snakemake@output[["pdf"]])[1] # Single PDF file path

# Parameters
PADJ_THRESH <- as.numeric(snakemake@params[["padj"]])[1] %||% 0.05
LFC_THRESH <- as.numeric(snakemake@params[["lfc"]])[1] %||% 0
TOP_N <- as.integer(snakemake@params[["top_n"]])[1] %||% 10

# ==========================================================================
# Processing & Plotting
# ==========================================================================
all_results_dt <- list()

message("INFO: Starting Volcano Plot generation for ", length(DESEQ_FILES), " files.")

# Open PDF device
pdf(PDF, width = 10, height = 8)

for (f in DESEQ_FILES) {
  if (!file.exists(f)) {
    message("WARNING: File does not exist: ", f)
    next
  }
  
  message("Processing file: ", basename(f))
  dds_obj <- readRDS(f)
  
  # Normalize input structure
  if (is.list(dds_obj) && !inherits(dds_obj, "DESeqDataSet") && !inherits(dds_obj, "DESeqResults")) {
    items_to_process <- dds_obj
  } else {
    items_to_process <- list(default = dds_obj)
  }

  for (analysis_name in names(items_to_process)) {
    dds <- items_to_process[[analysis_name]]
    if (is.null(dds)) next

    # CRITICAL FALLBACK: Detect structure type
    if (inherits(dds, "DESeqDataSet")) {
      contrasts <- resultsNames(dds)
      contrasts <- contrasts[contrasts != "Intercept"]
      message("  -> Found DESeqDataSet with ", length(contrasts), " contrasts.")
    } else {
      # If it is already a DESeqResults or a data.frame, create a fake single contrast loop
      contrasts <- "computed_results"
      message("  -> Found pre-computed results object. Processing directly.")
    }

    for (c in contrasts) {
      # Extract results table based on class
      if (inherits(dds, "DESeqDataSet")) {
        res_raw <- results(dds, name = c)
      } else {
        res_raw <- dds
      }
      
      # Convert to data.frame safely
      res_df <- as.data.frame(res_raw)
      
      # Check if required columns exist, otherwise skip to prevent empty plots
      if (!"padj" %in% colnames(res_df) || !"log2FoldChange" %in% colnames(res_df)) {
        message("  ERROR: Missing 'padj' or 'log2FoldChange' columns in object. Skipping.")
        next
      }
      
      # Convert to data.table
      res <- as.data.table(res_df, keep.rownames = "KO_Number")
      res[, `:=`(Analysis = analysis_name, Contrast = c)]
      
      # Set status tags
      res[, Color_Status := "Non Significatif"]
      res[padj <= PADJ_THRESH & log2FoldChange >= LFC_THRESH, Color_Status := "Sur"]
      res[padj <= PADJ_THRESH & log2FoldChange <= -LFC_THRESH, Color_Status := "Sous"]
      
      # Top labels
      res[, Label := ""]
      setorder(res, padj, na.last = TRUE)
      sig_rows <- which(!is.na(res$padj) & res$padj <= PADJ_THRESH)
      
      if (length(sig_rows) > 0) {
        top_rows <- head(sig_rows, TOP_N)
        res[top_rows, Label := KO_Number]
      }
      
      # Save to global list
      all_results_dt[[paste0(analysis_name, "_", c)]] <- res

      message("  -> Generating plot for contrast: ", c, " (Rows: ", nrow(res), ")")

      # Plotting
      p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Color_Status)) +
        geom_point(alpha = 0.4, size = 1.2) +
        geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(PADJ_THRESH), linetype = "dashed", alpha = 0.5) +
        geom_text_repel(aes(label = Label), size = 3, fontface = "bold", max.overlaps = 15) +
        scale_color_manual(
          values = c("Sur" = "forestgreen", "Sous" = "firebrick3", "Non Significatif" = "black"),
          drop = FALSE
        ) +
        theme_minimal() +
        labs(
          title = paste("Analysis:", analysis_name),
          subtitle = paste("Contrast:", c, "| padj <=", PADJ_THRESH),
          x = "Log2 Fold Change",
          y = "-log10(adj. P-value)"
        )
      
      print(p)
    }
  }
}

dev.off()

# ==========================================================================
# Final Data Export
# ==========================================================================
if (length(all_results_dt) > 0) {
  final_dt <- rbindlist(all_results_dt, use.names = TRUE, fill = TRUE)
  write_parquet(final_dt, PARQUET)
  message("✓ PDF report successfully updated: ", PDF)
} else {
  message("CRITICAL: all_results_dt is empty. No plots were drawn.")
}