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
library(ggrastr)
library(arrow)
library(rlang)
library(glue)

# Libraries Bioconductor
library(DESeq2)
library(phyloseq)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Inputs
DESEQ_FILES <- as.character(snakemake@input[["deseq_files"]])
PHYLOSEQ_OBJ <- as.character(snakemake@input[["phyloseq_obj"]])

# Outputs
PARQUET <- as.character(snakemake@output[["parquet"]])[1] # Single parquet file path
PDF <- as.character(snakemake@output[["pdf"]])[1] # Single PDF file path

# Shared plots features
SHARED <- snakemake@params[["shared"]]
THEME <- as.character(SHARED$theme) %||% "theme_minimal"
PDF_SIZE    <- as.numeric(SHARED$pdf_size) %||% c(14, 12)
TITLE_SIZE  <- as.integer(SHARED$title_size) %||% 12
SUBTITLE_SIZE <- as.integer(SHARED$subtitle_size) %||% 10
LEGEND_SIZE <- as.integer(SHARED$legend_size) %||% 10
AXES_SIZE   <- as.integer(SHARED$axes_size) %||% 10
RANK <- as.character(snakemake@params[["rank"]])[1]

# Parameters from config/params
TITLE_TEMPLATE <- as.character(snakemake@params[["title"]])[1] %||% "{source} | Top {top_n} Differential Abundance Volcano Plot in {contrast}"
SUBTITLE_TEMPLATE <- as.character(snakemake@params[["subtitle"]])[1] %||% "Pvalue threshold : {padj_threshold} | Log2 Fold Change threshold : {lfc_threshold} | {rank}"
CONTRAST_LIST <- unlist(as.character(snakemake@params[["contrast"]]))
PADJ_THRESHOLD <- as.numeric(unlist(snakemake@params[["padj"]])) %||% c(0.05, 0.05, 0.05)
LFC_THRESHOLD  <- as.numeric(unlist(snakemake@params[["lfc"]]))  %||% c(1.3, 1.3, 1.3)
TOP_N <- as.integer(snakemake@params[["top_n"]])[1] %||% 10

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ==========================================================================
# 1. Annotation Loading & Detection
# ==========================================================================
dt_annot <- NULL
is_kegg <- grepl("kegg", SOURCE, ignore.case = TRUE)

if (!is.na(PHYLOSEQ_OBJ) && file.exists(PHYLOSEQ_OBJ)) {
  message("INFO: Loading Phyloseq object for annotation: ", basename(PHYLOSEQ_OBJ))
  ps <- readRDS(PHYLOSEQ_OBJ)
  
  # On détermine dynamiquement le nom de l'identifiant selon la source
  target_id <- if (is_kegg) "kegg_id" else "contig_id"
  
  # Extraction en forçant les rownames de la tax_table à prendre le nom de notre target_id
  dt_annot <- as.data.table(as.data.frame(tax_table(ps)), keep.rownames = target_id)
}

# Si c'est du KEGG et que RANK n'est pas défini, on cible par défaut le niveau 3 (Pathways)
if (is_kegg) {
  KEGG_DISPLAY_RANK <- if (RANK != "") RANK else "level_3"
}

# ==========================================================================
# Processing & Plotting
# ==========================================================================
pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])

all_results_dt <- list()

for (i in seq_along(CONTRAST_LIST)) {
  
  contrast_type <- CONTRAST_LIST[i]
  CURRENT_PADJ <- PADJ_THRESHOLD[i]
  CURRENT_LFC  <- LFC_THRESHOLD[i]
  
  message(sprintf("\nProcessing contrast: %s", contrast_type))
  message(sprintf("  → Applying cutoffs: padj <= %g, lfc >= %g", CURRENT_PADJ, CURRENT_LFC))

  for (f in DESEQ_FILES) {
    if (!file.exists(f)) {
      message("WARNING: File does not exist: ", f)
      next
    }
    
    message("Processing model file: ", basename(f))
    master_rds <- readRDS(f) # Struct: list(models = list(...), dt = list(...))

    if (!contrast_type %in% names(master_rds$dt)) {
      message(sprintf("  ❌ [SKIP] Contrast type '%s' not found in $dt. (Available: %s)", 
                      contrast_type, paste(names(master_rds$dt), collapse = ", ")))
      next
    }
    
    res_full <- master_rds$dt[[contrast_type]]
    analysis_name <- paste0(tools::file_path_sans_ext(basename(f)), "_", contrast_type)
    
    if (is.null(res_full) || nrow(res_full) == 0) {
      message(sprintf("  ⚠️ [SKIP] Table for '%s' is empty.", contrast_type))
      next
    }
    
    # Extraire la liste de toutes les comparaisons générées pour ce modèle
    comparisons_to_plot <- unique(res_full$Comparison)
    message(sprintf("  → Found %d comparisons to plot for '%s'", length(comparisons_to_plot), contrast_type))
    
    # Loop through each comparison in the slot
    for (comp in comparisons_to_plot) {
      res <- res_full[Comparison == comp]
      
      if (!"padj" %in% names(res) || !"log2FoldChange" %in% names(res)) {
        next
      }
      
      target_id <- if (is_kegg) "kegg_id" else "contig_id"

      setnames(res, "Feature_ID", target_id, skip_absent = TRUE)

      # Jointure des annotations Phyloseq
      if (!is.null(dt_annot)) {
        res <- merge(res, dt_annot, by = target_id, all.x = TRUE)
        
        if (is_kegg) {
          # Traitement spécifique KEGG : Construction de l'étiquette combinée unique (ex: "K00163 | Pyruvate dehydrogenase")
          if (KEGG_DISPLAY_RANK %in% names(res)) {
            res[, Display_Name := {
              # Nettoyage des valeurs de description vides ou manquantes
              desc <- get(KEGG_DISPLAY_RANK)
              desc[is.na(desc) | desc == "" | desc == "Unassigned"] <- "Unknown Function"
              
              # On limite la taille de la description à 25 caractères pour ne pas surcharger le graphique
              desc_short <- substr(desc, 1, 25)
              
              # On crée le label combiné
              paste(get(target_id), desc_short, sep = " | ")
            }]
          } else {
            res[, Display_Name := get(target_id)]
          }
        } else {
          # Logique classique pour les contigs
          if (RANK %in% names(res)) {
            res[, Display_Name := get(RANK)]
          } else {
            res[, Display_Name := get(target_id)]
          }
        }
      } else {
        res[, Display_Name := get(target_id)]
      }

      # Remplacement des valeurs vides pour la sécurité
      res[is.na(Display_Name) | Display_Name == "", Display_Name := get(target_id)]
      
      # On limite la taille globale finale du label à 35 caractères max
      res[, Display_Name := substr(Display_Name, 1, 35)]
      
      # Assignation des couleurs (Sur / Sous / Non Significatif)
      res[, Color_Status := "Non Significatif"]
      res[padj <= CURRENT_PADJ & log2FoldChange >= CURRENT_LFC, Color_Status := "Sur"]
      res[padj <= CURRENT_PADJ & log2FoldChange <= -CURRENT_LFC, Color_Status := "Sous"]
      
      res[, Label := NA_character_]

      # Tri par p-value ajustée (les plus significatifs en premier)
      setorder(res, padj, na.last = TRUE)

      # On extrait temporairement les lignes significatives pour compter
      signif_rows <- which(res$padj <= CURRENT_PADJ)
      n_signif <- length(signif_rows)

      # On n'applique les labels QUE s'il y a au moins 1 gène significatif
      if (n_signif > 0) {
        # On détermine combien de gènes on va labelliser (au max TOP_N)
        n_to_label <- min(n_signif, TOP_N)
        
        # On cible les 'n_to_label' premières lignes significatives
        rows_to_label <- signif_rows[1:n_to_label]
        
        # Attribution du Display_Name à la colonne Label
        res[rows_to_label, Label := Display_Name]
      } else {
        message("⚠️ [Volcano] Aucun gène significatif trouvé (padj <= ", CURRENT_PADJ, "). Passage à la suite sans labels.")
      }
      
      message("  -> Drawing Volcano Plot for: ", comp)
      
      resolved_title <- glue(TITLE_TEMPLATE, source = toupper(SOURCE), top_n = TOP_N, contrast = toupper(comp))
      resolved_subtitle <- glue(SUBTITLE_TEMPLATE, padj_threshold = CURRENT_PADJ, lfc_threshold = CURRENT_LFC, rank = if (is_kegg) KEGG_DISPLAY_RANK else RANK)
      
      all_results_dt[[paste0(analysis_name, "_", comp)]] <- copy(res)
      
      theme_function <- match.fun(THEME)
      p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Color_Status)) +
        ggrastr::geom_point_rast(data = res[Color_Status == "Non Significatif"], alpha = 0.4, size = 1.2, raster.dpi = 150) +
        geom_point(data = res[Color_Status != "Non Significatif"], alpha = 0.4, size = 1.2) +
        geom_vline(xintercept = c(-CURRENT_LFC, CURRENT_LFC), linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(CURRENT_PADJ), linetype = "dashed", alpha = 0.5) +
        geom_text_repel(aes(label = Label), size = 3, fontface = "bold", max.overlaps = 15) +
        scale_color_manual(values = c("Sur" = "forestgreen", "Sous" = "firebrick3", "Non Significatif" = "black"), drop = FALSE) +
        labs(title = resolved_title, subtitle = resolved_subtitle, x = "Log2 Fold Change", y = "-log10(adj. P-value)") +
        theme_function() +
        theme(
          legend.position = "right",
          plot.title = element_text(size = TITLE_SIZE, face = "bold"),
          plot.subtitle = element_text(size = SUBTITLE_SIZE),
          legend.text = element_text(size = LEGEND_SIZE),
          axis.title = element_text(size = AXES_SIZE),
          axis.text = element_text(size = AXES_SIZE)
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