# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Volcano plot from Kegg data DESeq2
# PURPOSE : Differential Abundance Volcano Plots for DESeq2 Results
# AUTHOR  : Yann Le Bihan
# DATE    : 2026-09-03
# LINK    : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. METADATA & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
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
})

source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_DESEQ    <- as.character(snakemake@input[["deseq_files"]])[1]
IN_PHYLOSEQ <- as.character(snakemake@input[["phyloseq_obj"]])[1]

# Outputs
OUT_PARQUET <- as.character(snakemake@output[["parquet"]])[1]
OUT_PDF     <- as.character(snakemake@output[["pdf"]])[1]

# Shared Plot Parameters (Filtrés pour Volcano Plot)
PARAM_SHARED            <- snakemake@params[["shared"]]
PARAM_THEME             <- as.character(PARAM_SHARED$theme) %||% "theme_minimal"
PARAM_PDF_SIZE          <- as.numeric(PARAM_SHARED$pdf_size) %||% c(14, 12)
PARAM_TITLE_SIZE        <- as.numeric(PARAM_SHARED$title_size) %||% 12
PARAM_SUBTITLE_SIZE     <- as.numeric(PARAM_SHARED$subtitle_size) %||% 10
PARAM_AXES_TITLE_SIZE   <- as.numeric(PARAM_SHARED$axes_title_size) %||% 10
PARAM_AXES_TICK_SIZE    <- as.numeric(PARAM_SHARED$axes_tick_size) %||% 9
PARAM_LEGEND_TITLE_SIZE <- as.numeric(PARAM_SHARED$legend_title_size) %||% 10
PARAM_LEGEND_SIZE       <- as.numeric(PARAM_SHARED$legend_size) %||% 9
PARAM_POINT_SIZE        <- as.numeric(PARAM_SHARED$point_size) %||% 1.2
PARAM_VOLCANO_LABEL_SIZE<- as.numeric(PARAM_SHARED$volcano_label_size) %||% 3

# Specific Volcano Parameters & Translation Templates
TEMPLATE_TITLE    <- as.character(snakemake@params[["title_template"]])[1]
TEMPLATE_SUBTITLE <- as.character(snakemake@params[["subtitle_template"]])[1]

PARAM_RANK        <- tolower(as.character(snakemake@params[["rank"]])[1])
PARAM_CONTRAST    <- unlist(as.character(snakemake@params[["contrast"]]))
PARAM_PADJ        <- as.numeric(unlist(snakemake@params[["padj"]])) %||% c(0.05, 0.05, 0.05)
PARAM_LFC         <- as.numeric(unlist(snakemake@params[["lfc"]]))  %||% c(1.3, 1.3, 1.3)
PARAM_TOP_N       <- as.integer(snakemake@params[["top_n"]])[1] %||% 10

# Wildcards & Dynamic Variables
WILDCARD_SOURCE   <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ==========================================================================
# Helper : Résolution des textes/titres
# ==========================================================================
resolve_text <- function(template, vars = list()) {
  out <- template
  for (name in names(vars)) {
    out <- gsub(paste0("\\{", name, "\\}"), as.character(vars[[name]]), out)
  }
  return(out)
}

# Dynamic feature_type resolution
FEATURE_TYPE <- if (grepl("level_3", PARAM_RANK, ignore.case = TRUE)) {
  "pathways"
} else if (grepl("gene_description", PARAM_RANK, ignore.case = TRUE) || WILDCARD_SOURCE == "kegg") {
  "genes"
} else {
  "taxa"
}

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
# Validate the required files
for (f in c(IN_DESEQ, IN_PHYLOSEQ)) {
  if (!file.exists(f)) stop(sprintf("❌ Critical Error: Input file does not exist: %s", f))
}

# Verify the presence of RANK
if (is.na(PARAM_RANK) || PARAM_RANK == "" || PARAM_RANK == "null") {
  stop("❌ Critical Error: 'rank' parameter must be configured in config.yaml.")
}

# Source Validation
is_kegg   <- grepl("kegg", WILDCARD_SOURCE, ignore.case = TRUE)
is_contig <- grepl("contig", WILDCARD_SOURCE, ignore.case = TRUE)
is_reads  <- grepl("read", WILDCARD_SOURCE, ignore.case = TRUE)

if (!is_kegg && !is_contig && !is_reads) {
  stop(sprintf("❌ Critical Error: Unrecognized source '%s'. Must be kegg, contig, or reads.", WILDCARD_SOURCE))
}

TARGET_ID <- if (is_kegg) "kegg_id" else if (is_contig) "contig_id" else "read_id"

valid_contrasts <- c("ref", "combo", "date")

if (length(PARAM_CONTRAST) == 0 || 
    any(is.na(PARAM_CONTRAST)) || 
    !all(PARAM_CONTRAST %in% valid_contrasts)) {
  
  stop(sprintf(
    "❌ Critical Error: Invalid contrast parameter(s) [%s]. All values must be among: %s",
    paste(PARAM_CONTRAST, collapse = ", "),
    paste(valid_contrasts, collapse = ", ")
  ))
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading Phyloseq object...")
ps <- readRDS(IN_PHYLOSEQ)

if (!inherits(ps, "phyloseq")) stop("❌ Critical Error: Loaded object is not of class 'phyloseq'.")
if (phyloseq::nsamples(ps) == 0 || phyloseq::ntaxa(ps) == 0) stop("❌ Phyloseq object is empty.")

n_ps_samples <- phyloseq::nsamples(ps)
n_ps_taxa    <- phyloseq::ntaxa(ps)

if (is.null(n_ps_samples) || n_ps_samples == 0) {
  stop(sprintf("❌ Critical Error: Phyloseq object in '%s' contains 0 samples.", IN_PHYLOSEQ))
}

if (is.null(n_ps_taxa) || n_ps_taxa == 0) {
  stop(sprintf("❌ Critical Error: Phyloseq object in '%s' contains 0 taxa.", IN_PHYLOSEQ))
}

message(sprintf("✓ Phyloseq object loaded successfully: %d sample(s) and %d feature(s).", n_ps_samples, n_ps_taxa))

master_rds <- readRDS(IN_DESEQ)

if (is.null(master_rds) || !is.list(master_rds) || is.null(master_rds$dt) || length(master_rds$dt) == 0) {
  stop(sprintf("❌ Critical Error: DESeq2 object in '%s' is empty or invalid.", IN_DESEQ))
}
message(sprintf("✓ DESeq2 object loaded successfully: %d contrast(s) found (%s).", 
                length(master_rds$dt), paste(names(master_rds$dt), collapse = ", ")))

# Extracting Annotations for Joining
dt_annot <- if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
  dt <- as.data.table(as.data.frame(phyloseq::tax_table(ps)), keep.rownames = TARGET_ID)
  dt[, (TARGET_ID) := trimws(as.character(get(TARGET_ID)))]
  dt
} else {
  NULL
}

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PROCESSING LOOP
# ------------------------------------------------------------------------------
all_results_dt <- list()

with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  page_count <- 0

  for (i in seq_along(PARAM_CONTRAST)) {
    
    contrast_type <- PARAM_CONTRAST[i]
    current_padj  <- PARAM_PADJ[i]
    current_lfc   <- PARAM_LFC[i]
    
    message(sprintf("\nProcessing contrast: %s", contrast_type))
    message(sprintf("  → Applying cutoffs: padj <= %g, lfc >= %g", current_padj, current_lfc))
      
    if (!contrast_type %in% names(master_rds$dt)) {
      message(sprintf("  ❌ [SKIP] Contrast type '%s' not found in $dt. (Available: %s)", 
                      contrast_type, paste(names(master_rds$dt), collapse = ", ")))
      next
    }
    
    res_full <- master_rds$dt[[contrast_type]]
    analysis_name <- paste0(tools::file_path_sans_ext(basename(IN_DESEQ)), "_", contrast_type)
    
    if (is.null(res_full) || nrow(res_full) == 0) {
      message(sprintf("  ⚠️ [SKIP] Table for '%s' is empty.", contrast_type))
      next
    }
    
    comparisons_to_plot <- unique(res_full$Comparison)
    message(sprintf("  → Found %d comparisons to plot for '%s'", length(comparisons_to_plot), contrast_type))
    
    for (comp in comparisons_to_plot) {
      res <- copy(res_full[Comparison == comp])
      
      if (!"padj" %in% names(res) || !"log2FoldChange" %in% names(res)) {
        next
      }
      
      setnames(res, "Feature_ID", TARGET_ID, skip_absent = TRUE)

      # Joining / Formatting Annotations
      if (!is.null(dt_annot)) {
        res[, (TARGET_ID) := trimws(as.character(get(TARGET_ID)))]

        # Jointure avec dt_annot si PARAM_RANK n'est pas encore dans res
        if (!PARAM_RANK %in% names(res) && PARAM_RANK %in% names(dt_annot)) {
          res <- dt_annot[res, on = TARGET_ID]
        }

        # Définition de Display_Name
        if (PARAM_RANK %in% names(res)) {
          if (is_kegg) {
            # On travaille directement sur la colonne data.table
            res[, tmp_rank := get(PARAM_RANK)]
            res[is.na(tmp_rank) | tmp_rank == "" | tmp_rank == "Unassigned", tmp_rank := "Unknown Function"]
            res[, Display_Name := paste(get(TARGET_ID), substr(tmp_rank, 1, 25), sep = " | ")]
            res[, tmp_rank := NULL] # Nettoyage de la colonne temporaire
          } else {
            # Affichage classique du rang taxonomique
            res[, Display_Name := get(PARAM_RANK)]
          }
        } else {
          res[, Display_Name := get(TARGET_ID)]
        }
      } else {
        res[, Display_Name := get(TARGET_ID)]
      }
      
      res[is.na(Display_Name) | Display_Name == "", Display_Name := get(TARGET_ID)]
      res[, Display_Name := stringr::str_trunc(Display_Name, width = 35, side = "right")]
      
      # Assigning Color Status
      res[, Color_Status := "Non Significatif"]
      res[!is.na(padj) & padj <= current_padj & log2FoldChange >= current_lfc, Color_Status := "Sur"]
      res[!is.na(padj) & padj <= current_padj & log2FoldChange <= -current_lfc, Color_Status := "Sous"]
      
      res[, Label := NA_character_]
      setorder(res, padj, na.last = TRUE)
      
      signif_rows <- which(!is.na(res$padj) & res$padj <= current_padj)
      n_signif   <- length(signif_rows)
      actual_top_n <- min(n_signif, PARAM_TOP_N)
      
      if (n_signif > 0) {
        rows_to_label <- signif_rows[1:actual_top_n]
        res[rows_to_label, Label := Display_Name]
      } else {
        message("⚠️ [Volcano] Aucun gène significatif trouvé (padj <= ", current_padj, "). Passage à la suite sans labels.")
      }
      
      # ------------------------------------------------------------------------
      # 6. GRAPHICS GENERATION
      # ------------------------------------------------------------------------
      message("  -> Drawing Volcano Plot for: ", comp)
      
      resolved_title <- resolve_text(TEMPLATE_TITLE, list(
        source = toupper(WILDCARD_SOURCE),
        top_n = actual_top_n,
        feature_type = FEATURE_TYPE,
        contrast = toupper(comp)
      ))

      resolved_subtitle <- resolve_text(TEMPLATE_SUBTITLE, list(
        padj_threshold = current_padj,
        lfc_threshold = current_lfc,
        rank = PARAM_RANK
      ))
      
      all_results_dt[[paste0(analysis_name, "_", comp)]] <- copy(res)
      
      theme_function <- match.fun(PARAM_THEME)
      
      p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Color_Status)) +
        ggrastr::geom_point_rast(data = res[Color_Status == "Non Significatif"], alpha = 0.4, size = PARAM_POINT_SIZE, raster.dpi = 150) +
        geom_point(data = res[Color_Status != "Non Significatif"], alpha = 0.4, size = PARAM_POINT_SIZE) +
        geom_vline(xintercept = c(-current_lfc, current_lfc), linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(current_padj), linetype = "dashed", alpha = 0.5) +
        geom_text_repel(aes(label = Label), size = PARAM_VOLCANO_LABEL_SIZE, fontface = "bold", max.overlaps = 15) +
        scale_color_manual(values = c("Sur" = "forestgreen", "Sous" = "firebrick3", "Non Significatif" = "black"), drop = FALSE) +
        labs(title = resolved_title, subtitle = resolved_subtitle, x = "Log2 Fold Change", y = "-log10(adj. P-value)") +
        theme_function() +
        theme(
          legend.position  = "right",
          plot.title       = element_text(size = PARAM_TITLE_SIZE, face = "bold"),
          plot.subtitle    = element_text(size = PARAM_SUBTITLE_SIZE),
          axis.title       = element_text(size = PARAM_AXES_TITLE_SIZE),
          axis.text        = element_text(size = PARAM_AXES_TICK_SIZE),
          legend.title     = element_text(size = PARAM_LEGEND_TITLE_SIZE),
          legend.text      = element_text(size = PARAM_LEGEND_SIZE)
        )
      
      render_page(p)
      page_count <- page_count + 1
    }
  }

  if (page_count == 0) {
    render_fallback("No significant features found for Volcano plots.")
  }
})

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (length(all_results_dt) > 0) {
  final_dt <- rbindlist(all_results_dt, use.names = TRUE, fill = TRUE)
  arrow::write_parquet(final_dt, OUT_PARQUET)
  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  arrow::write_parquet(data.table(), OUT_PARQUET)
  message("⚠️ 'all_results_dt' is empty. No plots were drawn.")
}