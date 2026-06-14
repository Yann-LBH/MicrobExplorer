################################################################################
# Project : "MicrobExplorer"
# Script: "Utilitaires : Barplot_DESeq2"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
library(data.table)
library(ggplot2)
library(viridis)

# Libraries Bioconductor
library(phyloseq)

# ==========================================================================
# 1. Charger les résultats DESeq2
# ==========================================================================
load_deseq2_results <- function(deseq_paths) {
  #' Load DESeq2 results from RDS files
  #'
  #' @param deseq_paths Named list or named vector with file paths
  #'        Names should be contrast types (e.g., "ref", "date", "combo")
  #'        Values should be file paths
  #' @return List of data.tables with DESeq2 results

  results <- list()

  # Handle both named list and named vector
  if (!is.list(deseq_paths)) {
    deseq_paths <- as.list(deseq_paths)
  }

  for (contrast_type in names(deseq_paths)) {
    file_path <- deseq_paths[[contrast_type]]

    # Skip if file is NA or NULL
    if (is.na(file_path) || is.null(file_path)) {
      warning(sprintf("Skipping %s: path is NULL or NA", contrast_type))
      next
    }

    if (file.exists(file_path)) {
      tryCatch(
        {
          res <- readRDS(file_path)

          # Convert to data.table if it's a DESeq2Results object
          if (class(res)[1] == "DESeqResults") {
            res_dt <- as.data.table(as.data.frame(res), keep.rownames = "Feature_ID")
          } else if (is.data.frame(res) || is.data.table(res)) {
            res_dt <- as.data.table(res)
            # Check if Feature_ID exists, otherwise try to use rownames if available
            if (!"Feature_ID" %in% names(res_dt) && !"KO_Number" %in% names(res_dt)) {
              setnames(res_dt, 1, "Feature_ID")
            }
          } else {
            stop(sprintf("Unexpected class for %s: %s", contrast_type, class(res)[1]))
          }

          results[[contrast_type]] <- res_dt
          message(sprintf("✓ Loaded %s contrast: %d features", contrast_type, nrow(res_dt)))
        },
        error = function(e) {
          warning(sprintf("Failed to load %s: %s", contrast_type, e$message))
        }
      )
    } else {
      warning(sprintf("File not found for %s: %s", contrast_type, file_path))
    }
  }

  if (length(results) == 0) {
    stop("No DESeq2 results could be loaded")
  }

  return(results)
}

# ==========================================================================
# 2. Transformer phyloseq en data.table pour relative abundance
# ==========================================================================
phyloseq_to_dt <- function(ps_object) {
  #' Convert phyloseq object to data.table format
  #'
  #' @param ps_object phyloseq object (should already be relative abundance)
  #' @return data.table with columns: OTU, Sample, Abundance, Metadata...

  # Extract OTU table
  otu_dt <- as.data.table(as(otu_table(ps_object), "matrix"), keep.rownames = "OTU")

  # Extract taxonomy
  tax_dt <- as.data.table(as(tax_table(ps_object), "matrix"), keep.rownames = "OTU")

  # Extract metadata
  meta_dt <- as.data.table(as(sample_data(ps_object), "data.frame"), keep.rownames = "Sample")

  # Melt OTU table
  df <- melt(otu_dt, id.vars = "OTU", variable.name = "Sample", value.name = "Abundance")

  # Merge with taxonomy and metadata
  df <- merge(df, tax_dt, by = "OTU")
  df <- merge(df, meta_dt, by = "Sample")

  return(df)
}

# ==========================================================================
# 3. Obtenir les features significatives
# ==========================================================================
get_significant_features <- function(DESEQ_RESULTS, PADJ_THRESHOLD, LFC_THRESHOLD, direction) {
  #' Get significant features from DESeq2 results
  #'
  #' @param DESEQ_RESULTS data.table with DESeq2 results
  #' @param PADJ_THRESHOLD p-adjusted threshold (default 0.05)
  #' @param LFC_THRESHOLD log2 fold change threshold (default 0)
  #' @param direction "up" or "down" for over/under-represented
  #' @return Vector of feature IDs

  if (direction == "up") {
    features <- DESEQ_RESULTS[padj <= PADJ_THRESHOLD & log2FoldChange >= LFC_THRESHOLD, Feature_ID]
  } else if (direction == "down") {
    features <- DESEQ_RESULTS[padj <= PADJ_THRESHOLD & log2FoldChange <= -LFC_THRESHOLD, Feature_ID]
  } else {
    stop("direction must be 'up' or 'down'")
  }

  return(features)
}

# ==========================================================================
# 4. Générer la palette de couleurs
# ==========================================================================
generate_color_palette <- function(top_features, n_colors = TOP_N) {
  #' Generate color palette for top features
  #'
  #' @param top_features Vector of top feature names
  #' @param n_colors Number of distinct colors (default 10)
  #' @return Named vector of colors

  # Get color palette
  colors <- viridis::turbo(length(top_features))
  names(colors) <- top_features

  return(colors)
}

# ==========================================================================
# 5. Préparer les données pour le plot
# ==========================================================================
prepare_plot_data <- function(df_abundance, sig_features, TOP_N) {
  #' Prepare data for stacked barplot
  #'
  #' @param df_abundance data.table with abundance data
  #' @param sig_features Vector of significant feature IDs
  #' @param TOP_N Number of top features to highlight
  #' @return data.table formatted for plotting

  # Filter for significant features
  df_sig <- df_abundance[OTU %in% sig_features]

  # Get top features by abundance
  top_features <- df_sig[, .(Total_Abundance = sum(Abundance)), by = OTU][
    order(-Total_Abundance)
  ][1:min(TOP_N, .N), OTU]

  # Create category column
  df_abundance[, Category := "Other"]
  df_abundance[OTU %in% sig_features, Category := "Significant (Other)"]
  df_abundance[OTU %in% top_features, Category := OTU]

  # Aggregate by category and sample
  df_plot <- df_abundance[, .(Abundance = sum(Abundance)),
    by = .(Sample, Category)
  ]

  # Set factor levels
  factor_levels <- c("Other", "Significant (Other)", sort(top_features))
  df_plot[, Category := factor(Category, levels = factor_levels)]

  return(df_plot)
}

# ==========================================================================
# 6. Créer le plot stackedbarplot
# ==========================================================================
create_stackedbarplot <- function(df_plot, title, subtitle, group_column = "Sample") {
  #' Create stacked barplot
  #'
  #' @param df_plot data.table with columns: Sample, Category, Abundance
  #' @param title Plot title
  #' @param subtitle Plot subtitle
  #' @param group_column Column to group by (default "Sample")
  #' @return ggplot object

  # Generate colors
  categories <- unique(df_plot$Category)
  n_top <- sum(!categories %in% c("Other", "Significant (Other)"))

  color_map <- c(
    setNames(
      viridis::turbo(TOP_N),
      categories[!categories %in% c("Other", "Significant (Other)")]
    ),
    "Significant (Other)" = "grey60",
    "Other" = "grey90"
  )

  # Create plot
  p <- ggplot(df_plot, aes_string(x = group_column, y = "Abundance", fill = "Category")) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    scale_fill_manual(values = color_map, name = "Feature") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Sample",
      y = "Relative Abundance"
    )

  return(p)
}

# ==========================================================================
# Export functions
# ==========================================================================
export_results <- function(df_plot, output_path) {
  #' Export plot data as parquet
  #'
  #' @param df_plot data.table to export
  #' @param output_path Path to output file

  tryCatch(
    {
      arrow::write_parquet(df_plot, output_path)
      message(sprintf("✓ Results exported to: %s", output_path))
    },
    error = function(e) {
      warning(sprintf("Failed to export results: %s", e$message))
    }
  )
}
