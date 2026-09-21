# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Heatmap_from_DESeq2.R
# PURPOSE : Dynamic Heatmap for KEGG (KO & Pathways)
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
  library(rlang)
  library(circlize)
  library(viridis)
  library(vegan)
  library(arrow)
  # Libraries Bioconductor
  library(DESeq2)
  library(phyloseq)
  library(ComplexHeatmap)
})

source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
IN_DESEQ    <- as.character(snakemake@input[["deseq_files"]])[1]
IN_PHYLOSEQ <- as.character(snakemake@input[["phyloseq_obj"]])[1]

# Outputs
OUT_PDF     <- as.character(snakemake@output[["pdf"]])[1]
OUT_PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared Plot Parameters
PARAM_SHARED                <- snakemake@params[["shared"]]
PARAM_PALETTE               <- as.character(PARAM_SHARED$palette) %||% "turbo"
PARAM_PDF_SIZE              <- as.numeric(PARAM_SHARED$pdf_size) %||% c(7, 5.5)
PARAM_TITLE_SIZE            <- as.numeric(PARAM_SHARED$title_size) %||% 12
PARAM_LEGEND_TITLE_SIZE     <- as.numeric(PARAM_SHARED$legend_title_size) %||% 8.5
PARAM_LEGEND_SIZE           <- as.numeric(PARAM_SHARED$legend_size) %||% 8
PARAM_HEATMAP_LABEL_SIZE    <- as.numeric(PARAM_SHARED$heatmap_label_size) %||% 6.5
PARAM_DENDROGRAM_LINE_WIDTH <- as.numeric(PARAM_SHARED$dendrogram_line_width) %||% 0.6
PARAM_RANK                  <- as.character(snakemake@params[["rank"]])[1]

# Specific Heatmap DESeq2 Parameters
PARAM_TITLE             <- as.character(snakemake@params[["title_template"]])[1] %||% "{source} | Abundance of the top {top_n} taxons in sample {sample}"
PARAM_SUBTITLE          <- as.character(snakemake@params[["subtitle_template"]])[1] %||% "Hierarchical clustering: {clust_method} | Distance metric: {distance_method} | {rank}"
PARAM_CONTRAST          <- tolower(as.character(snakemake@params[["contrast"]]))
PARAM_PADJ              <- as.numeric(unlist(snakemake@params[["padj"]])) %||% c(0.05, 0.05, 0.05)
PARAM_LFC               <- as.numeric(unlist(snakemake@params[["lfc"]]))  %||% c(1.3, 1.3, 1.3)
PARAM_TOP_N             <- as.integer(snakemake@params[["top_n"]])[1] %||% 25
PARAM_GROUP_BY          <- as.character(snakemake@params[["group_by"]])[1] %||% "name"
PARAM_CLUST_METHOD      <- as.character(snakemake@params[["clust_method"]])[1] %||% "complete"
PARAM_DISTANCE_METHOD   <- as.character(snakemake@params[["distance_method"]])[1] %||% "bray"

# Wildcards & Variables globales
WILDCARD_SOURCE       <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

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
if (is.null(IN_DESEQ) || !file.exists(IN_DESEQ)) {
  stop(sprintf("❌ Critical Error: DESeq2 input RDS file '%s' does not exist.", IN_DESEQ))
}

if (is.null(IN_PHYLOSEQ) || !file.exists(IN_PHYLOSEQ)) {
  stop(sprintf("❌ Critical Error: Phyloseq input RDS file '%s' does not exist.", IN_PHYLOSEQ))
}

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
# 4. DATA LOADING
# ------------------------------------------------------------------------------
message("INFO: Loading input Deseq and Phyloseq...")
master_deseq <- readRDS(IN_DESEQ)
ps_obj       <- readRDS(IN_PHYLOSEQ)

# Extract metadata from phyloseq object
meta_df <- as.data.frame(sample_data(ps_obj))

if (!PARAM_GROUP_BY %in% colnames(meta_df)) {
  warning(sprintf("⚠️ WARNING: group_col '%s' not found in metadata. Falling back to 'name'.", PARAM_GROUP_BY))
  PARAM_GROUP_BY <- if ("name" %in% colnames(meta_df)) "name" else colnames(meta_df)[1]
}

KEGG_PATHWAY_RANKS <- c("level_1", "level_2", "level_3")
is_pathway_mode    <- any(sapply(KEGG_PATHWAY_RANKS, function(x) grepl(x, PARAM_RANK, ignore.case = TRUE)))

# Extraction nettoyée de la table de taxonomie
dt_tax_raw <- as.data.table(as.data.frame(tax_table(ps_obj)), keep.rownames = "Feature_ID")
setnames(dt_tax_raw, make.unique(colnames(dt_tax_raw)))
dt_tax_raw[, Feature_ID := trimws(as.character(Feature_ID))]

if (is_pathway_mode) {
  # Mode Pathway : on conserve l'annotation par voie métabolique
  message("INFO: Pathway mode detected. Preparing pathway taxonomy annotation...")
  dt_tax <- dt_tax_raw[, lapply(.SD, function(x) paste(unique(na.omit(x[x != ""])), collapse = "; ")), 
                       by = Feature_ID]
} else {
  # Mode KO / Gène : Si kegg_ids existe (cas d'un phyloseq composite), on déplie au besoin
  if ("kegg_ids" %in% names(dt_tax_raw)) {
    message("INFO: Unnesting kegg_ids to match KO-level DESeq2 features...")
    dt_tax <- dt_tax_raw[, .(kegg_id = trimws(unlist(strsplit(kegg_ids, "[;|\\,]")))), by = setdiff(names(dt_tax_raw), "kegg_ids")]
    setnames(dt_tax, "kegg_id", "Feature_ID")
  } else {
    dt_tax <- dt_tax_raw
  }
}

dt_tax[, Feature_ID := as.character(Feature_ID)]

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATION
# ------------------------------------------------------------------------------
message(sprintf("INFO: Processing DESeq2 results for contrast type: '%s'", PARAM_CONTRAST))

heatmap_data_list <- list()
filtered_dt_list  <- list()

for (i in seq_along(PARAM_CONTRAST)) {
    
  contrast_type <- PARAM_CONTRAST[i]
  current_padj  <- PARAM_PADJ[i]
  current_lfc   <- PARAM_LFC[i]

  dt_results <- master_deseq$dt[[contrast_type]]
  dds_model  <- master_deseq$models[[contrast_type]]

  if (is.null(dt_results) || nrow(dt_results) == 0 || is.null(dds_model)) {
    warning(sprintf("⚠️ No data found for contrast '%s', skipping...", contrast_type))
    next
  }

  dt_results[, Feature_ID := trimws(as.character(Feature_ID))]

  if (!PARAM_RANK %in% names(dt_results)) {
    n_before <- nrow(dt_results)
    dt_results <- merge(dt_results, dt_tax, by = "Feature_ID", all.x = TRUE)
    n_matched <- sum(!is.na(dt_results[[PARAM_RANK]]))
    message(sprintf("  → Taxonomy joined: %d/%d matched (%.1f%%)",
                    n_matched, n_before, 100 * n_matched / n_before))
    if (n_matched == 0) {
      warning(sprintf("⚠️ WARNING: 0 taxonomy matches for rank '%s' — check phyloseq/DESeq2 KO consistency.", PARAM_RANK))
    }
  }

  # Local helper for extracting matrices
  extract_heatmap_data <- function(dt, norm_counts, top_n, exclude_unassigned = FALSE) {
    taxa_col <- if (PARAM_RANK %in% names(dt)) PARAM_RANK else "Feature_ID"  
    work_dt <- data.table::copy(dt)

    if (exclude_unassigned) {
      work_dt <- work_dt[!get(taxa_col) %in% c("KO_Unassigned", "Pathway_Unassigned", "Unassigned", "") & 
                         !Feature_ID %in% c("KO_Unassigned", "Pathway_Unassigned", "Unassigned", "")]
    }
    
    sig_dt <- work_dt[!is.na(padj) & padj < current_padj & abs(log2FoldChange) >= current_lfc]
    sig_dt <- sig_dt[order(-abs(log2FoldChange))]
    top_dt <- head(sig_dt, top_n)
    
    if (nrow(top_dt) == 0) return(NULL)
    
    # Preserving Uniqueness by Feature_ID
    top_dt <- unique(top_dt, by = "Feature_ID")
    feat_ids <- top_dt$Feature_ID

    valid_feats <- intersect(feat_ids, rownames(norm_counts))
    if (length(valid_feats) == 0) {
      message("  ⚠️ Warning: Impossible d'aligner Feature_ID avec les lignes de norm_counts.")
      return(NULL)
    }
    mat_counts <- norm_counts[valid_feats, , drop = FALSE]
    top_dt     <- top_dt[Feature_ID %in% valid_feats]
    
    if (nrow(top_dt) == 0 || nrow(mat_counts) == 0) return(NULL)
    
    # Construction dynamique des étiquettes selon le mode
    if (is_pathway_mode) {
      rank_vals <- as.character(top_dt[[taxa_col]])
      rank_vals[is.na(rank_vals) | rank_vals == "" | rank_vals == "Unassigned"] <- "Unknown"
      # Pour un pathway : Affiche le nom du level_3
      top_dt[, feature_label := paste0(rank_vals, " (", Feature_ID, ")")]
    } else {
      # Pour un KO : Formatage de la description des gènes (Premier gène + (+ n ...))
      format_gene_desc <- function(raw_val) {
        if (is.na(raw_val) || raw_val == "" || raw_val == "Unassigned") return("Unknown")
        
        # 1. Découpage des blocs gènes séparés par une VIRGULE
        gene_blocks <- unlist(strsplit(as.character(raw_val), ","))
        
        # 2. Nettoyage de chaque bloc : suppression de ce qui suit le POINT-VIRGULE
        clean_genes <- sapply(gene_blocks, function(block) {
          trimws(sub(";.*$", "", block))
        })
        
        # Conservation des noms non vides
        clean_genes <- clean_genes[clean_genes != ""]
        
        # 3. Assemblage du label et décompte des gènes restants
        if (length(clean_genes) > 1) {
          return(sprintf("%s (+ %d ...)", clean_genes[1], length(clean_genes) - 1))
        } else if (length(clean_genes) == 1) {
          return(clean_genes[1])
        } else {
          return(raw_val)
        }
      }

      top_dt[, formatted_rank := sapply(get(taxa_col), format_gene_desc)]

      top_dt[, feature_label := fifelse(
        Feature_ID == formatted_rank,
        Feature_ID,
        paste(Feature_ID, formatted_rank, sep = " | ")
      )]
      
      top_dt[, formatted_rank := NULL]
    }
    
    scaled_mat <- t(scale(t(mat_counts)))
    scaled_mat[is.na(scaled_mat)] <- 0
    return(list(scaled_mat = scaled_mat, top_dt = top_dt))
  }

  norm_counts <- counts(dds_model, normalized = TRUE)
    
  for (comp in unique(dt_results$Comparison)) {
    comp_dt <- dt_results[Comparison == comp]
      
    data_with_unassigned <- extract_heatmap_data(comp_dt, norm_counts, PARAM_TOP_N, exclude_unassigned = FALSE)
    data_no_unassigned   <- extract_heatmap_data(comp_dt, norm_counts, PARAM_TOP_N, exclude_unassigned = TRUE)
      
    # Unique key that prevents overlap between two contrasts
    key_name <- paste0(contrast_type, "_", comp)

    if (!is.null(data_with_unassigned) || !is.null(data_no_unassigned)) {
      heatmap_data_list[[key_name]] <- list(
        with_unassigned = data_with_unassigned,
        no_unassigned   = data_no_unassigned
      )
        
      if (!is.null(data_with_unassigned)) filtered_dt_list[[paste0(comp, "_with")]] <- data_with_unassigned$top_dt
      if (!is.null(data_no_unassigned))   filtered_dt_list[[paste0(comp, "_no")]]   <- data_no_unassigned$top_dt
    }
  }
}

# ------------------------------------------------------------------------------
# 6. VISUALIZATION / OUTPUT
# ------------------------------------------------------------------------------
if (length(heatmap_data_list) == 0) {
  warning(sprintf("⚠️ WARNING: No significant features (padj < %.2f) found for contrast type '%s'.", 
                  PARAM_PADJ, PARAM_CONTRAST))
}

# Helper function to render a single heatmap
draw_heatmap_instance <- function(mat_data, top_df, comp_title, subtitle_str, meta_df) {
  scaled_max  <- max(abs(mat_data), na.rm = TRUE)
  if (scaled_max == 0) scaled_max <- 1
  col_fun_mat <- colorRamp2(c(-scaled_max, 0, scaled_max), c("blue", "white", "red"))

  top_df_std <- as.data.frame(top_df)
  idx <- match(as.character(rownames(mat_data)), as.character(top_df_std$Feature_ID))
  top_df_aligned <- top_df_std[idx, , drop = FALSE]
  lfc_vec        <- top_df_aligned$log2FoldChange

  row_labels_vec <- top_df_aligned$feature_label
  if (is.null(row_labels_vec)) row_labels_vec <- rownames(mat_data)

  lfc_max     <- max(abs(lfc_vec), na.rm = TRUE)
  if (lfc_max == 0) lfc_max <- 1
  col_fun_lfc <- colorRamp2(c(-lfc_max, 0, lfc_max), c("darkgreen", "white", "darkorange"))
  
  sub_meta <- meta_df[colnames(mat_data), , drop = FALSE]
  anno_col_name <- if (PARAM_GROUP_BY %in% colnames(meta_df)) {
    PARAM_GROUP_BY
  } else {
    warning(sprintf(
      "⚠️ WARNING: group_col '%s' not found in metadata (%s). Falling back to first available column.",
      PARAM_GROUP_BY, paste(colnames(meta_df), collapse = ", ")
    ))
    colnames(meta_df)[1]
  }
  group_vals    <- as.character(sub_meta[[anno_col_name]])

  # Annotations avec tailles de police et titres appliqués
  col_anno <- HeatmapAnnotation(
    Group = group_vals,
    col   = list(Group = setNames(viridis(length(unique(group_vals))), 
                                unique(group_vals))),
    annotation_legend_param = list(
      Group = list(
        title_gp = grid::gpar(fontsize = PARAM_LEGEND_TITLE_SIZE, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = PARAM_LEGEND_SIZE)
      )
    )
  )
  
  row_anno <- rowAnnotation(
    Log2FC = lfc_vec,
    col    = list(Log2FC = col_fun_lfc),
    annotation_legend_param = list(
      Log2FC = list(
        title_gp = grid::gpar(fontsize = PARAM_LEGEND_TITLE_SIZE, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = PARAM_LEGEND_SIZE)
      )
    )
  )
  full_title <- paste0(comp_title, "\n", subtitle_str)
  ht <- Heatmap(
    matrix              = mat_data,
    name                = "Z-Score",
    col                 = col_fun_mat,
    top_annotation      = col_anno,
    left_annotation     = row_anno,
    show_row_names      = TRUE,
    row_labels          = top_df_aligned$feature_label,
    row_names_gp        = grid::gpar(fontsize = PARAM_HEATMAP_LABEL_SIZE, fontitalic = TRUE),
    row_names_max_width = grid::unit(12, "cm"),
    show_column_names   = TRUE,
    column_names_gp     = grid::gpar(fontsize = PARAM_HEATMAP_LABEL_SIZE),
    cluster_rows        = TRUE,
    cluster_columns     = TRUE,
    # Style des dendrogrammes (épaisseur de ligne)
    row_dend_gp         = grid::gpar(lwd = PARAM_DENDROGRAM_LINE_WIDTH),
    column_dend_gp      = grid::gpar(lwd = PARAM_DENDROGRAM_LINE_WIDTH),
    # Titre principal et sous-titre stylisés
    column_title        = full_title,
    column_title_gp     = grid::gpar(fontsize = PARAM_TITLE_SIZE, fontface = "bold"),
    # Style de la légende principale (Z-Score)
    heatmap_legend_param = list(
      title_gp  = grid::gpar(fontsize = PARAM_LEGEND_TITLE_SIZE, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = PARAM_LEGEND_SIZE)
    )
  )
  
  render_page(ht)
}

is_pathway <- exists("PARAM_RANK") && grepl("pathway|level", PARAM_RANK, ignore.case = TRUE)
unassigned_label <- if (is_pathway) "Pathway_Unassigned" else "KO_Unassigned"

with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  page_count <- 0
  
  if (length(heatmap_data_list) > 0) {
    for (comp in names(heatmap_data_list)) {
      comp_entry <- heatmap_data_list[[comp]]
      
      # 1. Plot avec éléments non assignés
      if (!is.null(comp_entry$with_unassigned)) {
        title_str <- resolve_text(PARAM_TITLE, list(
          source = toupper(WILDCARD_SOURCE),
          top_n = nrow(comp_entry$with_unassigned$top_dt),
          feature_type = if (toupper(WILDCARD_SOURCE) == "KEGG") "genes" else "taxa",
          comparison = comp,
          unassigned_status = sprintf("(With %s)", unassigned_label)
        ))
        subtitle_str <- resolve_text(PARAM_SUBTITLE, list(
          contrast        = comp, 
          clust_method    = PARAM_CLUST_METHOD,
          distance_method = PARAM_DISTANCE_METHOD,
          padj_threshold  = PARAM_PADJ[1],
          lfc_threshold   = PARAM_LFC[1],
          rank            = PARAM_RANK
        ))
        draw_heatmap_instance(comp_entry$with_unassigned$scaled_mat, comp_entry$with_unassigned$top_dt, title_str, subtitle_str, meta_df)
        page_count <- page_count + 1
      }
      
      # 2. Plot sans éléments non assignés
      if (!is.null(comp_entry$no_unassigned)) {
        title_str <- resolve_text(PARAM_TITLE, list(
          source = toupper(WILDCARD_SOURCE),
          top_n = nrow(comp_entry$no_unassigned$top_dt),
          feature_type = if (toupper(WILDCARD_SOURCE) == "KEGG") "genes" else "taxa",
          comparison = comp,
          unassigned_status = sprintf("(Excl. %s)", unassigned_label)
        ))
        subtitle_str <- resolve_text(PARAM_SUBTITLE, list(
          contrast        = comp,
          clust_method    = PARAM_CLUST_METHOD,
          distance_method = PARAM_DISTANCE_METHOD,
          padj_threshold  = PARAM_PADJ[1],
          lfc_threshold   = PARAM_LFC[1],
          rank            = PARAM_RANK
        ))
        draw_heatmap_instance(comp_entry$no_unassigned$scaled_mat, comp_entry$no_unassigned$top_dt, title_str, subtitle_str, meta_df)
        page_count <- page_count + 1
      }
    }
  }

  if (page_count == 0) {
    render_fallback("No significant features found for Heatmaps.")
  }
})

# ------------------------------------------------------------------------------
# 7. EXPORT & CLEANUP
# ------------------------------------------------------------------------------
if (length(filtered_dt_list) > 0) {
  master_parquet <- rbindlist(filtered_dt_list, use.names = TRUE, fill = TRUE)
  arrow::write_parquet(master_parquet, OUT_PARQUET)
  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Parquet : ", OUT_PARQUET)
} else {
# Write empty parquet structure if no data processed
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: Empty output: No results written.")
}