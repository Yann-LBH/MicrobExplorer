# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Heatmap.R
# PURPOSE : Dynamic Heatmap for Taxonomy (Reads & Contigs)
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
  library(phyloseq)
  library(ComplexHeatmap)
})

source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_PHYLOSEQ <- as.character(snakemake@input[["phyloseq_obj"]])[1]

# Outputs
OUT_PDF     <- as.character(snakemake@output[["pdf"]])[1]
OUT_PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared Plot Parameters
PARAM_SHARED                <- snakemake@params[["shared"]]
PARAM_PALETTE               <- as.character(PARAM_SHARED$palette) %||% "turbo"
PARAM_PDF_SIZE              <- as.numeric(PARAM_SHARED$pdf_size) %||% c(7, 5.5)
PARAM_TITLE_SIZE            <- as.numeric(PARAM_SHARED$title_size) %||% 12
PARAM_AXES_TITLE_SIZE       <- as.numeric(PARAM_SHARED$axes_title_size) %||% 9.5
PARAM_AXES_TICK_SIZE        <- as.numeric(PARAM_SHARED$axes_tick_size) %||% 8
PARAM_LEGEND_TITLE_SIZE     <- as.numeric(PARAM_SHARED$legend_title_size) %||% 8.5
PARAM_LEGEND_SIZE           <- as.numeric(PARAM_SHARED$legend_size) %||% 8
PARAM_HEATMAP_LABEL_SIZE    <- as.numeric(PARAM_SHARED$heatmap_label_size) %||% 6.5
PARAM_DENDROGRAM_LINE_WIDTH <- as.numeric(PARAM_SHARED$dendrogram_line_width) %||% 0.6
PARAM_STROKE_WIDTH          <- as.numeric(PARAM_SHARED$stroke_width) %||% 0.5

PARAM_RANK                  <- as.character(snakemake@params[["rank"]])[1]

# Specific Heatmap Parameters & Translation Templates
TEMPLATE_TITLE    <- as.character(snakemake@params[["title_template"]])[1]
TEMPLATE_SUBTITLE <- as.character(snakemake@params[["subtitle_template"]])[1]
TEXT_SAMPLE_ALL   <- as.character(snakemake@params[["text_sample_all"]])[1]

PARAM_TOP_N             <- as.integer(snakemake@params[["top_n"]])[1] %||% 50
PARAM_CLUST_METHOD      <- as.character(snakemake@params[["clust_method"]])[1] %||% "complete"
PARAM_DISTANCE_METHOD   <- as.character(snakemake@params[["distance_method"]])[1] %||% "bray"

# Wildcards & Variables globales
WILDCARD_SOURCE         <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# ==========================================================================
# Helper : Resolution des textes/titres
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

# Helper for dynamic transposition and calculation tool based on the distance method
compute_taxa_distance <- function(mat, method = "bray") {
  method <- tolower(method)
  
  # Ecological metrics using vegan::vegdist (requires taxa in columns -> t(mat))
  vegan_methods <- c("bray", "canberra", "kulczynski", "jaccard", "gower", "altgower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis")
  
  # Standard metrics via stats::dist (requires taxa in lines -> mat)
  dist_methods <- c("euclidean", "maximum", "manhattan", "binary", "minkowski")

  if (method %in% vegan_methods) {
    d <- vegan::vegdist(mat, method = method)
  } else if (method %in% dist_methods) {
    d <- stats::dist(mat, method = method)
  } else {
    stop(sprintf("❌ Unsupported distance method: '%s'. Please select one of the following: %s", 
                 method, paste(c(vegan_methods, dist_methods), collapse = ", ")))
  }
  
  d_mat <- as.matrix(d)
  d_mat[is.na(d_mat) | is.nan(d_mat)] <- 0
  return(as.dist(d_mat))
}

# ==========================================================================
# 1. Chargement de l'objet Phyloseq
# ==========================================================================
ps <- readRDS(IN_PHYLOSEQ)

# ==========================================================================
# 2. Agglomération selon le RANK (Strict : Erreur si RANK invalide)
# ==========================================================================
if (PARAM_RANK != "") {
  available_ranks <- rank_names(ps)
  
  if (!PARAM_RANK %in% available_ranks) {
    stop(sprintf(
      "❌ Configuration Error: The requested RANK '%s' was not found in phyloseq object ranks [%s].",
      PARAM_RANK, paste(available_ranks, collapse = ", ")
    ))
  }
  
  message("--- Aggregating at rank: ", PARAM_RANK, " ---")
  ps <- tax_glom(ps, taxrank = PARAM_RANK, NArm = FALSE)
}

# ==========================================================================
# 3. Échelle de couleurs GLOBALE (Calculée après agglomération)
# ==========================================================================
global_mat <- as.matrix(otu_table(ps)@.Data)
if (!taxa_are_rows(ps)) global_mat <- t(global_mat)
global_mat[is.na(global_mat)] <- 0
global_max <- max(log10(global_mat + 1), na.rm = TRUE)
if (is.na(global_max) || global_max <= 0) global_max <- 1

col_fun <- colorRamp2(seq(0, global_max, length.out = 5), viridis(5, option = tolower(PARAM_PALETTE)))

# ==========================================================================
# 4. Génération Heatmaps (PDF) + Export Tidy Data (Shiny)
# ==========================================================================
samples <- unique(sample_data(ps)$name)
if (anyNA(samples)) {
  warning(sprintf("%d sample(s) with NA 'name' in metadata will be skipped.", sum(is.na(samples))))
  samples <- samples[!is.na(samples)]
}

# Helper function to generate a single heatmap
generate_heatmap <- function(ps_sub, s) {
  ps_work <- ps_sub
  
  # Identify target taxonomy column
  tax_mat_all <- as.matrix(tax_table(ps_work))
  target_rank <- if (PARAM_RANK != "") PARAM_RANK else colnames(tax_mat_all)[ncol(tax_mat_all)]
  
  if (ntaxa(ps_work) == 0) return(NULL)
  
  # Select Top N taxa
  top_taxa <- names(sort(taxa_sums(ps_work), decreasing = TRUE))[seq_len(min(PARAM_TOP_N, ntaxa(ps_work)))]
  ps_top   <- prune_taxa(top_taxa, ps_work)
  
  mat <- as.matrix(otu_table(ps_top)@.Data)
  if (!taxa_are_rows(ps_top)) mat <- t(mat)
  
  mat[is.na(mat)] <- 0
  mat <- log10(mat + 1)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(NULL)
  
  ps_top      <- prune_taxa(rownames(mat), ps_top)
  tax_mat_top <- as.matrix(tax_table(ps_top))
  
  if (!target_rank %in% colnames(tax_mat_top)) {
    stop(sprintf("❌ Configuration Error: Target rank '%s' not present in tax_table columns.", target_rank))
  }
  
  tax_vals <- as.character(tax_mat_top[, target_rank])
  tax_vals[is.na(tax_vals) | tax_vals == "" | tax_vals == "Unassigned"] <- "Unknown"

  # Drawing Heatmap
  cluster_rows_param <- FALSE
  if (nrow(mat) > 1) {
    cluster_rows_param <- tryCatch({
      dist_obj <- compute_taxa_distance(mat, method = PARAM_DISTANCE_METHOD)
      hclust(dist_obj, method = PARAM_CLUST_METHOD)
    }, error = function(e) {
      warning("⚠️ Distance computation/hclust failed. Row clustering disabled: ", e$message)
      return(FALSE)
    })
  }

  # Dynamic title and subtitle generation
  title_str <- resolve_text(TEMPLATE_TITLE, list(
    source = toupper(WILDCARD_SOURCE),
    top_n = nrow(mat),
    feature_type = FEATURE_TYPE,
    sample = if (!is.null(s) && s != "") s else TEXT_SAMPLE_ALL
  ))

  subtitle_str <- resolve_text(TEMPLATE_SUBTITLE, list(
    clust_method = PARAM_CLUST_METHOD,
    distance_method = PARAM_DISTANCE_METHOD,
    rank = PARAM_RANK
  ))

  full_title <- paste0(title_str, "\n", subtitle_str)

  ht <- Heatmap(
      mat,
      column_title = full_title,
      column_title_gp = grid::gpar(fontsize = PARAM_TITLE_SIZE, fontface = "bold"),
      name = "Abondance\n(log10)",
      heatmap_legend_param = list(
        title_gp = grid::gpar(fontsize = PARAM_LEGEND_TITLE_SIZE, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = PARAM_LEGEND_SIZE),
        legend_height = grid::unit(5, "cm"),
        grid_width = grid::unit(0.8, "cm")
      ),
      cluster_rows = cluster_rows_param,
      cluster_columns = FALSE,
      row_dend_width = grid::unit(25, "mm"),
      row_dend_gp = grid::gpar(lwd = PARAM_DENDROGRAM_LINE_WIDTH),
      row_labels = tax_vals,
      row_names_gp = grid::gpar(fontsize = PARAM_HEATMAP_LABEL_SIZE, fontitalic = TRUE),
      column_names_gp = grid::gpar(fontsize = PARAM_AXES_TICK_SIZE),
      col = col_fun
  )
  
  ht_drawn <- ComplexHeatmap::draw(ht)
  render_page(ht_drawn)

  # Extraction for Shiny
  row_idx          <- row_order(ht_drawn)
  mat_ordered      <- mat[row_idx, , drop = FALSE]
  tax_vals_ordered <- tax_vals[row_idx]

  dt_shiny  <- as.data.table(as.table(mat_ordered))
  setnames(dt_shiny, c("feature_id", "condition_col", "log10_abundance"))
  label_map <- setNames(tax_vals_ordered, rownames(mat_ordered))

  dt_shiny[, `:=`(
    sample_name   = s,
    feature_label = label_map[as.character(feature_id)],
    raw_abundance = (10^log10_abundance) - 1,
    rank          = target_rank,
    source        = WILDCARD_SOURCE,
    row_order     = rep(seq_along(row_idx), times = ncol(mat_ordered))
  )]
  
  return(dt_shiny)
}

# Dynamic PDF layout setup based on output volume
with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  page_count <- 0

  heatmap_loop <- rbindlist(lapply(samples, function(s) {
    keep_samples <- sample_names(ps)[sample_data(ps)$name == s]
    ps_sub       <- prune_samples(keep_samples, ps)
    
    if (ntaxa(ps_sub) == 0) {
      render_fallback(paste("Pas assez de données pour :", s))
      page_count <<- page_count + 1
      return(NULL)
    }
    
    # Generate single heatmap per sample
    dt_res <- generate_heatmap(ps_sub, s)
    
    if (!is.null(dt_res)) {
      page_count <<- page_count + 1
      return(dt_res)
    } else {
      render_fallback(paste("Pas assez de données pour :", s))
      page_count <<- page_count + 1
      return(NULL)
    }
  }), fill = TRUE)

  if (page_count == 0) {
    render_fallback("No features found for Heatmaps.")
  }
})

# ==========================================================================
# 5. Sauvegarde du fichier Parquet consolidé pour Shiny
# ==========================================================================

if (!is.null(heatmap_loop) && nrow(heatmap_loop) > 0) {
  arrow::write_parquet(heatmap_loop, OUT_PARQUET)
  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  # Write empty parquet structure if no data processed
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: Empty output: No results written.")
}