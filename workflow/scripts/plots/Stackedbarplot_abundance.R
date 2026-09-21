# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Stackedbarplot_abundance.R
# PURPOSE : Dynamic Stacked Barplots for Taxonomy & KEGG Pathways (Absolute/Relative)
# AUTHOR  : Yann Le Bihan
# DATE    : 2025-12-01
# LINK    : https://github.com/Yann-LBH/MicrobExplorer
#
# MODES (SNAKEMAKE$PARAMS$MODE):
#   "pathway"   — level_3 categories from KEGG pathway files
#   "taxonomy"  — taxonomic rank from intersec + contigs sources
#   "organisms" — genus-level RPKM, produces global bar + stacked bar
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. METADATA & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Libraries CRAN
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(rlang)
  library(arrow)
})

source("workflow/scripts/utils/utils_io.R")
source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_PHYLOSEQ     <- as.character(snakemake@input[["phyloseq_obj"]])[1]

# Outputs
OUT_PDF         <- as.character(snakemake@output[["pdf"]])[1]
OUT_PARQUET     <- as.character(snakemake@output[["parquet"]])[1]

# Shared Plot Parameters (Filtrés pour Stacked Barplot)
PARAM_SHARED            <- snakemake@params[["shared"]]
PARAM_THEME             <- as.character(PARAM_SHARED$theme) %||% "theme_minimal"
PARAM_PALETTE           <- as.character(PARAM_SHARED$palette) %||% "turbo"
PARAM_PDF_SIZE          <- c(10, 8)  # as.numeric(PARAM_SHARED$pdf_size) %||% c(12, 8)
PARAM_TITLE_SIZE        <- as.numeric(PARAM_SHARED$title_size) %||% 14
PARAM_SUBTITLE_SIZE     <- as.numeric(PARAM_SHARED$subtitle_size) %||% 10
PARAM_AXES_TITLE_SIZE   <- as.numeric(PARAM_SHARED$axes_title_size) %||% 10
PARAM_AXES_TICK_SIZE    <- as.numeric(PARAM_SHARED$axes_tick_size) %||% 9
PARAM_LEGEND_TITLE_SIZE <- as.numeric(PARAM_SHARED$legend_title_size) %||% 10
PARAM_LEGEND_SIZE       <- as.numeric(PARAM_SHARED$legend_size) %||% 9

# Specific Stackedbarplot Parameters & Translation Templates
TEMPLATE_TITLE    <- as.character(snakemake@params[["title_template"]])[1]
TEMPLATE_SUBTITLE <- as.character(snakemake@params[["subtitle_template"]])[1]
TEXT_SCOPE_ALL    <- as.character(snakemake@params[["text_scope_all"]])[1]
TEXT_SCOPE_EACH   <- as.character(snakemake@params[["text_scope_each"]])[1]

PARAM_MODE        <- as.character(snakemake@params[["mode"]])[1]
PARAM_TOP_N       <- as.integer(snakemake@params[["top_n"]])[1] %||% 10
PARAM_STAND_COL   <- tolower(as.character(snakemake@params[["stand_col"]])[1])
PARAM_RANK        <- as.character(snakemake@params[["rank"]])[1]

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

phyloseq_to_dt_fast <- function(ps, value_col = PARAM_STAND_COL) {
  # 1. Extraction et formatage de la matrice OTU/Abondance
  otu_mat <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }
  
  dt_long <- data.table::as.data.table(otu_mat, keep.rownames = "Feature_ID")
  dt_long <- data.table::melt(
    dt_long, 
    id.vars = "Feature_ID", 
    variable.name = "Sample", 
    value.name = value_col
  )
  
  # 2. Fusion avec les métadonnées échantillons (ex: name, date)
  if (!is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))) {
    dt_meta <- data.table::as.data.table(
      as(phyloseq::sample_data(ps), "data.frame"), 
      keep.rownames = "Sample"
    )
    dt_long <- merge(dt_long, dt_meta, by = "Sample", all.x = TRUE)
  }
  
  # 3. Fusion avec la table de taxonomie
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    dt_tax <- data.table::as.data.table(
      as(phyloseq::tax_table(ps), "matrix"), 
      keep.rownames = "Feature_ID"
    )
    dt_long <- merge(dt_long, dt_tax, by = "Feature_ID", all.x = TRUE)
  }
  
  return(dt_long)
}

color_palette <- function(categories, palette_name = "turbo") {
  top_cats  <- sort(setdiff(categories, "Others"))
  lvl_order <- c(top_cats, "Others")
  colours   <- c(viridisLite::viridis(length(top_cats), option = tolower(palette_name)), "#000000")
  names(colours) <- lvl_order
  return(list(colours = colours, levels = lvl_order))
}

stacked_bar <- function(dt, x_col, y_col, fill_col, colours, title, subtitle,
                        x_lab, y_lab, fill_lab, facet_col = NULL,
                        theme_name, title_size, 
                        subtitle_size, legend_size,
                        axes_title_size, axes_tick_size, 
                        legend_title_size, stroke_width){
  theme_function <- match.fun(theme_name)
  p <- ggplot(dt, aes(x = as.factor(.data[[x_col]]), y = .data[[y_col]], fill = .data[[fill_col]])) +
    geom_bar(
      stat = "identity",
      position = position_stack(reverse = TRUE)
    ) +
    scale_fill_manual(values = colours, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab, fill = fill_lab) +
    theme_function() +
    theme(
      plot.title         = element_text(size = title_size, face = "bold"),
      plot.subtitle      = element_text(size = subtitle_size, face = "italic"),
      axis.title         = element_text(size = axes_title_size),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = axes_tick_size),
      axis.text.y        = element_text(size = axes_tick_size),
      legend.title       = element_text(size = legend_title_size),
      legend.text        = element_text(size = legend_size - 2),
      panel.grid.major.x = element_blank()
    )
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)))
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
if (is.null(PARAM_STAND_COL) || is.na(PARAM_STAND_COL) || PARAM_STAND_COL == "") {
  stop("❌ Critical Error: 'stand_col' parameter is missing or empty in Snakemake config.")
}

for (f in IN_PHYLOSEQ) {
  if (!file.exists(f)) stop(sprintf("❌ Critical Error: Input file does not exist: %s", f))
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading Phyloseq object...")
ps <- readRDS(IN_PHYLOSEQ)
if (!inherits(ps, "phyloseq")) stop("❌ Critical Error: Loaded object is not of class 'phyloseq'.")

dt_merged <- phyloseq_to_dt_fast(ps, value_col = PARAM_STAND_COL)

# Cleaning up unassigned values in the extracted taxonomy
if (PARAM_RANK %in% names(dt_merged)) {
  dt_merged[is.na(get(PARAM_RANK)) | get(PARAM_RANK) == "", (PARAM_RANK) := "Unclassified"]
} else {
  stop(sprintf("❌ Critical Error: Rank column [%s] not found in Phyloseq object.", PARAM_RANK))
}

if (!"Feature_ID" %in% names(dt_merged)) {
  stop("❌ ERROR: 'Feature_ID' column missing from dt_merged — check phyloseq_to_dt_fast() output.")
}

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PROCESSING FUNCTIONS
# ------------------------------------------------------------------------------

# Mode 1: Relative Abundance of Pathways (KEGG)
run_pathway_relative_global <- function() {
  
  # Création dynamique de la colonne combinée (KO | Description) pour KEGG
  is_kegg <- isTRUE(get0("is_kegg")) || grepl("kegg", get0("WILDCARD_SOURCE"), ignore.case = TRUE)
  target_feature_col <- PARAM_RANK
  
  if (is_kegg) {
    target_feature_col <- "combined_label"
    
    found_col <- intersect("Feature_ID", names(dt_merged))
    target_ko_col <- if (length(found_col) > 0) found_col[1] else "taxa_id"
    
    # Nettoyage et sécurisation des valeurs
    ko_vals <- as.character(dt_merged[[target_ko_col]])
    ko_vals[is.na(ko_vals) | ko_vals == "" | ko_vals == "NA"] <- "KO_Unassigned"
    
    rank_vals <- as.character(dt_merged[[PARAM_RANK]])
    rank_vals[is.na(rank_vals) | rank_vals == "" | rank_vals == "Unassigned"] <- "Unclassified"
    
    # Fonction locale pour extraire le 1er nom et compter les autres dans la description
    format_kegg_label <- function(ko, desc) {
      if (ko == "KO_Unassigned" || desc == "Unclassified") return(paste(ko, desc, sep = " | "))
      
      # 1. Découpage des différents gènes séparés par une VIRGULE
      gene_blocks <- unlist(strsplit(desc, ","))
      
      # 2. Pour chaque bloc gène, on nettoie en supprimant ce qui suit le POINT-VIRGULE
      clean_genes <- sapply(gene_blocks, function(block) {
        trimws(sub(";.*$", "", block))
      })
      
      # Conservation des noms non vides
      clean_genes <- clean_genes[clean_genes != ""]
      
      # 3. Assemblage du label principal et décompte des autres gènes (+ n ...)
      if (length(clean_genes) > 1) {
        return(sprintf("%s | %s (+ %d ...)", ko, clean_genes[1], length(clean_genes) - 1))
      } else if (length(clean_genes) == 1) {
        return(sprintf("%s | %s", ko, clean_genes[1]))
      } else {
        return(sprintf("%s | %s", ko, desc))
      }
    }

    # Application de la mise en forme ligne par ligne (ou via Map)
    dt_merged[, (target_feature_col) := Map(format_kegg_label, ko_vals, rank_vals)]
    dt_merged[, (target_feature_col) := as.character(get(target_feature_col))]
  }

  if (!target_feature_col %in% names(dt_merged)) {
    stop(sprintf("The column [%s] is missing.", target_feature_col))
  }

  # Dynamic Subtitle Resolution
  subtitle_str <- resolve_text(TEMPLATE_SUBTITLE, list(
    mode = PARAM_MODE,
    stand_col = PARAM_STAND_COL,
    rank = PARAM_RANK
  ))

  # 1. Compute global abundance ONCE to define Top N (including Unassigned)
  global_agg <- dt_merged[, .(Global_Abund = sum(get(PARAM_STAND_COL), na.rm = TRUE)), by = c(target_feature_col)][order(-Global_Abund)]
  actual_top_n <- min(PARAM_TOP_N, nrow(global_agg))
  top_cats   <- global_agg[seq_len(actual_top_n)][[target_feature_col]]

  # Determine full category order with "Others" at the end
  all_possible_cats <- unique(c(top_cats, "Others"))
  
  # MASTER PALETTE: Generated once for all categories (including Unassigned)
  pal_master <- color_palette(all_possible_cats, palette_name = PARAM_PALETTE)

  # --------------------------------------------------------------------------
  # --- Plot 1: Global Horizontal Barplot (Includes Unassigned) ---
  # --------------------------------------------------------------------------
  title_all <- resolve_text(TEMPLATE_TITLE, list(
    source = toupper(WILDCARD_SOURCE),
    top_n = actual_top_n,
    feature_type = FEATURE_TYPE,
    scope = TEXT_SCOPE_ALL
  ))

  plot1_dt <- global_agg[, .(category = fifelse(get(target_feature_col) %in% top_cats, get(target_feature_col), "Others"), Global_Abund)][
    , .(Total_Abund = sum(Global_Abund)), by = category][order(Total_Abund)]
  
  cat_order_p1 <- plot1_dt[order(Total_Abund)]$category
  cat_order_p1 <- c(setdiff(cat_order_p1, "Others"), "Others")
  plot1_dt[, category := factor(category, levels = cat_order_p1)]

  theme_function <- match.fun(PARAM_THEME)

  p_global <- ggplot(plot1_dt, aes(x = Total_Abund, y = category, fill = category)) +
    geom_col(show.legend = FALSE) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
    theme_function() + 
    theme(
      plot.title   = element_text(face = "bold", size = PARAM_TITLE_SIZE),
      plot.subtitle= element_text(size = PARAM_SUBTITLE_SIZE),
      axis.title   = element_text(size = PARAM_AXES_TITLE_SIZE),
      axis.text.x  = element_text(size = PARAM_AXES_TICK_SIZE),
      axis.text.y  = element_text(size = PARAM_AXES_TICK_SIZE, face = "italic")
    ) +
    labs(
      title = title_all, 
      subtitle = subtitle_str, 
      x = paste("Total", toupper(PARAM_STAND_COL)), 
      y = PARAM_RANK
    ) +
    scale_fill_manual(values = pal_master$colours, drop = FALSE)

  # --------------------------------------------------------------------------
  # --- Plot 2: Relative Stacked Barplot (Excludes Unassigned, same palette) ---
  # --------------------------------------------------------------------------
  title_each <- resolve_text(TEMPLATE_TITLE, list(
    source = toupper(WILDCARD_SOURCE),
    top_n = actual_top_n,
    feature_type = FEATURE_TYPE,
    scope = TEXT_SCOPE_EACH
  ))

  # Filter out Unassigned entries
  dt_filtered <- dt_merged[!grepl("Unclassified", get(target_feature_col), ignore.case = TRUE)]

  plot2_dt <- dt_filtered[, .(total = sum(get(PARAM_STAND_COL), na.rm = TRUE)), 
                          by = c("name", "date", target_feature_col)]
  
  # Assign top categories based on the initial global Top N
  plot2_dt[, category := fifelse(get(target_feature_col) %in% top_cats, get(target_feature_col), "Others")]
  
  final_dt <- plot2_dt[, .(sum_val = sum(total)), by = .(name, date, category)][
    , pct := (sum_val / sum(sum_val)) * 100, by = .(name, date)]
  
  # Preserve initial factor ordering for visual consistency
  cats_p2_levels <- intersect(cat_order_p1, unique(final_dt$category))
  final_dt[, category := factor(category, levels = cats_p2_levels)]

  p_stacked <- stacked_bar(
    final_dt, "date", "pct", "category", pal_master$colours, title_each, subtitle_str, 
    "Date", "Relative Abundance (%)", "Pathways", facet_col = "name",
    theme_name = PARAM_THEME, title_size = PARAM_TITLE_SIZE, 
    subtitle_size = PARAM_SUBTITLE_SIZE, axes_title_size = PARAM_AXES_TITLE_SIZE,
    axes_tick_size = PARAM_AXES_TICK_SIZE, legend_title_size = PARAM_LEGEND_TITLE_SIZE,
    legend_size = PARAM_LEGEND_SIZE
  )

  render_page(p_global)
  render_page(p_stacked)
  return(2)
}

# Mode 2 : Relative Abundance per Sample
run_relative_by_sample <- function() {

  agg <- dt_merged[, .(Abund_Sum = sum(get(PARAM_STAND_COL), na.rm = TRUE)), by = c("name", "date", PARAM_RANK)]
  setnames(agg, PARAM_RANK, "Taxon")
  
  agg[, Abund_Pct := (Abund_Sum / sum(Abund_Sum)) * 100, by = .(name, date)]
  
  top_taxa_table <- agg[, .(G = sum(Abund_Sum)), by = Taxon][order(-G)]
  actual_top_n   <- min(PARAM_TOP_N, nrow(top_taxa_table))
  top_taxa       <- top_taxa_table[seq_len(actual_top_n), Taxon]

  agg[, Taxon_Final := fifelse(Taxon %in% top_taxa, Taxon, "Others")]

  final_dt  <- agg[, .(Abund_Pct = sum(Abund_Pct)), by = .(name, date, Taxon_Final)]
  tax_order <- final_dt[, .(total_abundance = sum(Abund_Pct)), by = Taxon_Final][order(total_abundance)]$Taxon_Final
  tax_order <- c(setdiff(tax_order, "Others"), "Others")
    
  pal <- color_palette(unique(final_dt$Taxon_Final), palette_name = PARAM_PALETTE)
  final_dt[, Taxon_Final := factor(Taxon_Final, levels = tax_order)]

  title_each <- resolve_text(TEMPLATE_TITLE, list(
    source = toupper(WILDCARD_SOURCE),
    top_n = actual_top_n,
    feature_type = FEATURE_TYPE,
    scope = TEXT_SCOPE_EACH
  ))

  subtitle_str <- resolve_text(TEMPLATE_SUBTITLE, list(
    mode = PARAM_MODE,
    stand_col = PARAM_STAND_COL,
    rank = PARAM_RANK
  ))

  p <- stacked_bar(final_dt, "date", "Abund_Pct", "Taxon_Final", pal$colours, title_each, subtitle_str, 
                   "Date", "Relative Abundance (%)", PARAM_RANK, facet_col = "name",
                   theme_name = PARAM_THEME, title_size = PARAM_TITLE_SIZE, 
                   subtitle_size = PARAM_SUBTITLE_SIZE, axes_title_size = PARAM_AXES_TITLE_SIZE,
                   axes_tick_size = PARAM_AXES_TICK_SIZE, legend_title_size = PARAM_LEGEND_TITLE_SIZE,
                   legend_size = PARAM_LEGEND_SIZE)
  render_page(p)
  return(1)
}

# Mode 3 : Total Absolute Abundance + Per Sample
run_absolute_global <- function() {
  if (!PARAM_RANK %in% names(dt_merged)) stop(sprintf("The taxonomy column [%s] is missing.", PARAM_RANK))
  dt_merged[get(PARAM_RANK) == "" | is.na(get(PARAM_RANK)), (PARAM_RANK) := "Unclassified"]

  taxa_config_global <- dt_merged[, .(Global_Abund = sum(get(PARAM_STAND_COL), na.rm = TRUE)), by = c(PARAM_RANK)][order(-Global_Abund)]
  actual_top_n       <- min(PARAM_TOP_N, nrow(taxa_config_global))
  top_genera         <- taxa_config_global[seq_len(actual_top_n)] [[PARAM_RANK]]

  subtitle_str <- resolve_text(TEMPLATE_SUBTITLE, list(
    mode = PARAM_MODE,
    stand_col = PARAM_STAND_COL,
    rank = PARAM_RANK
  ))

  # Plot 1: Barplot horizontal global
  title_all <- resolve_text(TEMPLATE_TITLE, list(
    source = toupper(WILDCARD_SOURCE),
    top_n = actual_top_n,
    feature_type = FEATURE_TYPE,
    scope = TEXT_SCOPE_ALL
  ))

  plot1_dt   <- taxa_config_global[, .(Taxa_Grouped = fifelse(get(PARAM_RANK) %in% top_genera, get(PARAM_RANK), "Others"), Global_Abund)][, .(Total_Abund = sum(Global_Abund)), by = Taxa_Grouped][order(Total_Abund)]
  taxa_order <- plot1_dt[order(Total_Abund)]$Taxa_Grouped
  taxa_order <- c(setdiff(taxa_order, "Others"), "Others")
  plot1_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = taxa_order)]

  pal_global     <- color_palette(levels(plot1_dt$Taxa_Grouped), palette_name = PARAM_PALETTE)
  theme_function <- match.fun(PARAM_THEME)
    
  p_global <- ggplot(plot1_dt, aes(x = Total_Abund, y = Taxa_Grouped, fill = Taxa_Grouped)) +
    geom_col(show.legend = FALSE) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
    theme_function() + 
    theme(
      plot.title   = element_text(face = "bold", size = PARAM_TITLE_SIZE),
      plot.subtitle= element_text(size = PARAM_SUBTITLE_SIZE),
      axis.title   = element_text(size = PARAM_AXES_TITLE_SIZE),
      axis.text.x  = element_text(size = PARAM_AXES_TICK_SIZE),
      axis.text.y  = element_text(size = PARAM_AXES_TICK_SIZE, face = "italic")
    ) +
    labs(title = title_all, subtitle = subtitle_str, x = paste("Total", toupper(PARAM_STAND_COL)), y = PARAM_RANK) +
    scale_fill_manual(values = pal_global$colours, drop = FALSE)

  # Plot 2: Stacked Barplot by Sample
  title_each <- resolve_text(TEMPLATE_TITLE, list(
    source = toupper(WILDCARD_SOURCE),
    top_n = actual_top_n,
    feature_type = FEATURE_TYPE,
    scope = TEXT_SCOPE_EACH
  ))

  plot2_dt <- dt_merged[!is.na(name) & !is.na(date), .(Total_Abund = sum(get(PARAM_STAND_COL), na.rm = TRUE)), by = c("name", "date", PARAM_RANK)]
  plot2_dt[, Taxa_Grouped := fifelse(get(PARAM_RANK) %in% top_genera, get(PARAM_RANK), "Others")]
  plot2_dt <- plot2_dt[, .(Total_Abund = sum(Total_Abund)), by = .(name, date, Taxa_Grouped)]
  plot2_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = taxa_order)]

  p_stacked <- stacked_bar(plot2_dt, "date", "Total_Abund", "Taxa_Grouped", pal_global$colours, title_each, subtitle_str, 
                           "Sampling Date", "Abundance", PARAM_RANK, facet_col = "name",
                           theme_name = PARAM_THEME, title_size = PARAM_TITLE_SIZE, 
                           subtitle_size = PARAM_SUBTITLE_SIZE, axes_title_size = PARAM_AXES_TITLE_SIZE,
                           axes_tick_size = PARAM_AXES_TICK_SIZE, legend_title_size = PARAM_LEGEND_TITLE_SIZE,
                           legend_size = PARAM_LEGEND_SIZE)

  render_page(p_global)
  render_page(p_stacked)
  return(2)
}

# ------------------------------------------------------------------------------
# 6. EXECUTION CORE & GRAPHICS GENERATION
# ------------------------------------------------------------------------------
message(sprintf("INFO: Executing Stacked Barplot module in mode '%s'...", PARAM_MODE))

with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  page_count <- switch(PARAM_MODE,
    "pathway_relative"   = run_pathway_relative_global(),
    "relative_by_sample"        = run_relative_by_sample(),
    "absolute_global"           = run_absolute_global(),
    stop("❌ Critical Error: Unknown mode specified: ", PARAM_MODE)
  )
  if (is.null(page_count) || page_count == 0) {
    render_fallback("⚠️ WARNING: Empty output: No results written..")
  }
})

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (length(dt_merged) && nrow(dt_merged) > 0) {
  arrow::write_parquet(dt_merged, OUT_PARQUET)
  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: Empty output: No results written.")
}