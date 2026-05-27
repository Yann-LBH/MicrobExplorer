################################################################################
# Project : "MicrobExplorer"
# Script: "Volcano plot from Kegg data DESeq2"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(ggplot2)
library(ggrepel)
library(arrow)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================
# Inputs from Snakemake
rds_files      <- snakemake@input[["rds_files"]]   # List of RDS files
parquet_output <- snakemake@output[["parquet"]]   # Single parquet file path
pdf_output     <- snakemake@output[["pdf"]]       # Single PDF file path

# Thresholds from rule params
padj_thresh <- snakemake@params[["padj"]] %||% 0.05
lfc_thresh  <- snakemake@params[["lfc"]]  %||% 0

# ==========================================================================
# Processing & Plotting
# ==========================================================================
all_results_dt <- list()

# Open PDF device to capture all subsequent plots
pdf(pdf_output, width = 10, height = 8)

for (f in rds_files) {
  dds <- readRDS(f)
  file_name <- basename(f)
  analysis_name <- sub("^(deseq2_)?(.*)\\.rds$", "\\2", file_name)
  contrast_label <- NA_character_
  if (grepl("_(ref|date|combo)(?:_.*)?\\.rds$", file_name)) {
    contrast_label <- sub(".*_(ref|date|combo)(?:_.*)?\\.rds$", "\\1", file_name)
  }
  
  process_result <- function(res, contrast) {
    res <- as.data.table(as.data.frame(res), keep.rownames = "KO_Number")
    res[, `:=`(Analysis = analysis_name, Contrast = contrast)]
    res[, Color_Status := "Non Significatif"]
    res[padj <= padj_thresh & log2FoldChange >= lfc_thresh, Color_Status := "Sur"]
    res[padj <= padj_thresh & log2FoldChange <= -lfc_thresh, Color_Status := "Sous"]
    res[, Label := ""]
    setorder(res, padj, na.last = TRUE)
    res[padj <= padj_thresh][1:10, Label := KO_Number]
    all_results_dt[[paste0(analysis_name, "_", contrast)]] <<- res
    
    p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Color_Status)) +
      geom_point(alpha = 0.4, size = 1.2) +
      geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", alpha = 0.5) +
      geom_text_repel(aes(label = Label), size = 3, fontface = "bold", max.overlaps = 15) +
      scale_color_manual(
        values = c("Sur" = "forestgreen", "Sous" = "firebrick3", "Non Significatif" = "black")
      ) +
      theme_minimal() +
      labs(
        title = paste("Analysis:", analysis_name),
        subtitle = paste("Contrast:", contrast, "| padj <=", padj_thresh),
        x = "Log2 Fold Change",
        y = "-log10(adj. P-value)"
      )
    print(p)
  }
  
  if (inherits(dds, "DESeqDataSet")) {
    contrasts <- resultsNames(dds)
    contrasts <- contrasts[contrasts != "Intercept"]
    for (c in contrasts) {
      res <- results(dds, name = c)
      process_result(res, c)
    }
  } else if (inherits(dds, "DESeqResults")) {
    process_result(dds, ifelse(is.na(contrast_label), "unknown", contrast_label))
  } else {
    stop(sprintf("Unsupported input type for file %s: %s", f, class(dds)[1]))
  }
}

dev.off() # Close PDF

# ==========================================================================
# Final Data Export
# ==========================================================================
# Combine all results into one table and save as Parquet
final_dt <- rbindlist(all_results_dt)
write_parquet(final_dt, parquet_output)

message("✓ PDF report generated: ", pdf_output)
message("✓ Data exported to Parquet: ", parquet_output)