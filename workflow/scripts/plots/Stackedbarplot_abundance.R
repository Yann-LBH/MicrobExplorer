################################################################################
# Project : "MicrobExplorer"
# Script  : "Unified stacked barplot — 
#            pathway relative / relative by sample / absolute global"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
#
# Modes (snakemake$params$mode):
#   "pathway"   — level_3 categories from KEGG pathway files
#   "taxonomy"  — taxonomic rank from intersec + contigs sources
#   "organisms" — genus-level RPKM, produces global bar + stacked bar
################################################################################

# ==========================================================================
# Snakemake configuration
# ==========================================================================

# Libraries CRAN
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(rlang)
  library(arrow)
})
# ==========================================================================
# Configuration Snakemake
# ==========================================================================
source("workflow/scripts/utils/utils_title_resolver.R")
source("workflow/scripts/utils/utils_io.R")


# Inputs
DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])[1]
TAXONOMY <- as.character(snakemake@input[["taxonomy"]])[1]

# Outputs
PDF     <- as.character(snakemake@output[["pdf"]])[1]
PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared plots features
SHARED      <- snakemake@params[["shared"]]
THEME       <- as.character(SHARED$theme) %||% "theme_minimal"
PALETTE     <- as.character(SHARED$palette) %||% "turbo"
PDF_SIZE    <- as.numeric(SHARED$pdf_size) %||% c(12, 8)
TITLE_SIZE  <- as.integer(SHARED$title_size) %||% 14
SUBTITLE_SIZE <- as.integer(SHARED$subtitle_size) %||% 10
LEGEND_SIZE <- as.integer(SHARED$legend_size) %||% 10
AXES_SIZE   <- as.integer(SHARED$axes_size) %||% 10

# Parameters
MODE              <- as.character(snakemake@params[["mode"]])[1]
TOP_N             <- as.integer(snakemake@params[["top_n"]])[1] %||% 10
STAND_COL         <- tolower(as.character(snakemake@params[["stand_col"]])[1])
RANK              <- as.character(snakemake@params[["rank"]])[1]

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# Standardization of Structure Type
FEATURE_TYPE <- if (SOURCE == "kegg") "pathways" else "taxon"

# TITRES ET SUBTITLES UNIFIÉS
RESOLVED_TITLE    <- list(
  all = resolve_text("TITLE_STACKEDBARPLOT", source = toupper(SOURCE), top_n = TOP_N, feature_type = FEATURE_TYPE, scope = t("all")),
  each = resolve_text("TITLE_STACKEDBARPLOT", source = toupper(SOURCE), top_n = TOP_N, feature_type = FEATURE_TYPE, scope = t("each"))
)
RESOLVED_SUBTITLE <- resolve_text("SUBTITLE_STACKEDBARPLOT", mode = MODE, stand_col = toupper(STAND_COL), rank = RANK)

# ==========================================================================
# Data loading
# ==========================================================================

# 1. Metadata loading
meta_dt <- load_metadata(METADATA)

id_col <- switch(SOURCE,
  "reads"   = "read_id",
  "contigs" = "contig_id",
  "kegg"    = "kegg_id"
)

# 2. Optimized TSV loading (MUST include id_col so left_join doesn't fail)
all_data <- load_tsv_dir_dynamic(
  paths       = DATA,
  meta_dt     = meta_dt,
  select_cols = unique(c(id_col, "sample_id", RANK, STAND_COL))
)

# ==========================================================================
# TAXONOMY loading and merging
# ==========================================================================
dt_taxo <- fread(TAXONOMY)

# 3. Conditional merge with deduplication using data.table syntax
if (SOURCE == "reads") {
  dt_taxo_clean <- unique(dt_taxo, by = "tax_id")
  dt_merged     <- merge(all_data, dt_taxo_clean, by.x = "read_id", by.y = "tax_id", all.x = TRUE)
  tax_ranks     <- c("tax_id", "scientific_name", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

} else if (SOURCE == "contigs") {
  dt_taxo_clean <- unique(dt_taxo, by = "contig_id")
  dt_merged     <- merge(all_data, dt_taxo_clean, by.x = "contig_id", by.y = "contig_id", all.x = TRUE)
  tax_ranks     <- c("domain", "phylum", "class", "order", "family", "genus", "species")

} else { # KEGG or other sources
  dt_taxo_clean <- unique(dt_taxo, by = "ko")
  dt_merged     <- merge(all_data, dt_taxo_clean, by.x = "kegg_id", by.y = "ko", all.x = TRUE)
  tax_ranks     <- c("ec_number", "level_1", "level_2", "level_3", "gene_description")
}

# 4. Clean up unassigned taxonomy in-place (replaces mutate + across)
tax_ranks_present <- intersect(tax_ranks, names(dt_merged))
for (col in tax_ranks_present) {
  set(dt_merged, i = which(is.na(dt_merged[[col]])), j = col, value = "Unclassified")
}
dt_merged <- dt_merged[meta_dt, on = "sample_id", nomatch = 0L]
# ==========================================================================
# Global Graphic initialization
# ==========================================================================

color_palette <- function(categories) {
  top_cats <- sort(setdiff(categories, "Others"))
  lvl_order <- c(top_cats, "Others")
  colours <- c(viridisLite::viridis(length(top_cats), option = tolower(PALETTE)), "#000000")
  names(colours) <- lvl_order
  list(colours = colours, levels = lvl_order)
}

# Universal stacked bar
stacked_bar <- function(dt, x_col, y_col, fill_col, colours, title, subtitle,
                        x_lab, y_lab, fill_lab, facet_col = NULL) {
  theme_function <- match.fun(THEME)
  p <- ggplot(dt, aes(x = as.factor(get(x_col)), y = get(y_col), fill = get(fill_col))) +
    geom_bar(
      stat = "identity",
      position = position_stack(reverse = TRUE),
      colour = "white",
      linewidth = 0.05
    ) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab, fill = fill_lab) +
    theme_function() +
    theme(
      plot.title = element_text(size = TITLE_SIZE, face = "bold"),
      plot.subtitle = element_text(size = SUBTITLE_SIZE, face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = LEGEND_SIZE - 2),
      panel.grid.major.x = element_blank()
    )
  
  if (is.list(colours)) colours <- colours$colours
  p <- p + scale_fill_manual(values = colours)
  
  # Dynamic Facet Addition (Used in run_absolute_global mode)
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_col)), scales = "free_x") +
      theme(strip.text = element_text(face = "bold", size = 12), panel.spacing = unit(1, "lines"))
  }
  return(p)
}

# ==========================================================================
# RUN MODES
# ==========================================================================

run_pathway_relative <- function(subtitle = RESOLVED_SUBTITLE) {
  target_col <- STAND_COL
  if (is.na(target_col) || target_col == "") stop("L'argument 'stand_col' est manquant ou égal à NA.")

  write_parquet(dt_merged, PARQUET)

  agg <- dt_merged[, .(total = sum(get(target_col), na.rm = TRUE)), by = c("name", "date", RANK)]
  agg[, category := fifelse(frank(-total, ties.method = "first") <= TOP_N, as.character(get(RANK)), "Others"), by = .(name, date)]
  
  final_dt <- agg[, .(sum_val = sum(total)), by = .(name, date, category)][, pct := (sum_val / sum(sum_val)) * 100, by = .(name, date)]
  # Stack Order Management
  cat_order <- final_dt[, .(total_abundance = sum(sum_val)), by = category][order(total_abundance)]$category
  cat_order <- c(setdiff(cat_order, "Others"), "Others") # keep other at top
  pal <- color_palette(unique(final_dt$category))
  final_dt[, category := factor(category, levels = cat_order)]

  pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
  on.exit(if (names(dev.cur()) != "null device") dev.off())

  # Using dynamic title matching the exact sample name
  print(stacked_bar(final_dt, "date", "pct", "category", pal, RESOLVED_TITLE$each, subtitle, "Date", "Relative Abundance (%)", "Pathways", facet_col = "name"))
}

run_relative_by_sample <- function(subtitle = RESOLVED_SUBTITLE) {
  target_col <- STAND_COL
  if (is.na(target_col) || target_col == "") stop("L'argument 'stand_col' est manquant ou égal à NA.")

  write_parquet(dt_merged, PARQUET)

  setnames(dt_merged, RANK, "Taxon")
  agg <- dt_merged[, .(Abund_Sum = sum(get(target_col), na.rm = TRUE)), by = c("name", "date", "Taxon")][, Abund_Pct := (Abund_Sum / sum(Abund_Sum)) * 100, by = .(name, date)]
  top_taxa <- agg[, .(G = sum(Abund_Sum)), by = Taxon][order(-G)[seq_len(min(TOP_N, .N))], Taxon]
  agg[, Taxon_Final := fifelse(Taxon %in% top_taxa, Taxon, "Others")]

  final_dt <- agg[, .(Abund_Pct = sum(Abund_Pct)), by = .(name, date, Taxon_Final)]
  # Stack Order Management
  tax_order <- final_dt[, .(total_abundance = sum(Abund_Pct)), by = Taxon_Final][order(total_abundance)]$Taxon_Final
  tax_order <- c(setdiff(tax_order, "Others"), "Others")
  pal <- color_palette(unique(final_dt$Taxon_Final))
  final_dt[, Taxon_Final := factor(Taxon_Final, levels = tax_order)]
  
  pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
  on.exit(if (names(dev.cur()) != "null device") dev.off(), add = TRUE)

  # Standardized title structure to match the each view
  print(stacked_bar(final_dt, "date", "Abund_Pct", "Taxon_Final", pal, RESOLVED_TITLE$each, subtitle, "Date", "Relative Abundance (%)", RANK, facet_col = "name"))
}

run_absolute_global <- function(subtitle = RESOLVED_SUBTITLE) {
  target_col <- STAND_COL
  if (is.na(target_col) || target_col == "") stop("L'argument 'stand_col' est manquant ou égal à NA.")

  if (!RANK %in% names(dt_merged)) stop(sprintf("The taxonomy column [%s] is missing.", RANK))
  dt_merged[get(RANK) == "" | is.na(get(RANK)), (RANK) := "Unclassified"]
  write_parquet(dt_merged, PARQUET)

  taxa_config_global <- dt_merged[, .(Global_Abund = sum(get(target_col), na.rm = TRUE)), by = c(RANK)][order(-Global_Abund)]
  top_genera <- taxa_config_global[seq_len(min(TOP_N, .N))][[RANK]]

  # --- Plot 1: Global horizontal bar ---
  plot1_dt <- taxa_config_global[, .(Taxa_Grouped = fifelse(get(RANK) %in% top_genera, get(RANK), "Others"), Global_Abund)][, .(Total_Abund = sum(Global_Abund)), by = Taxa_Grouped][order(Total_Abund)]
  # Stack Order Management
  taxa_order <- plot1_dt[order(Total_Abund)]$Taxa_Grouped
  taxa_order <- c(setdiff(taxa_order, "Others"), "Others")
  plot1_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = taxa_order)]

  # --- 3. Palette générée APRÈS la création de l'ordre ---
  pal_global <- color_palette(levels(plot1_dt$Taxa_Grouped))
  
  theme_function <- match.fun(THEME)
  p_global <- ggplot(plot1_dt, aes(x = Total_Abund, y = Taxa_Grouped, fill = Taxa_Grouped)) +
    geom_col(show.legend = FALSE) + 
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
    theme_function() + 
    theme(plot.title = element_text(face = "bold", size = TITLE_SIZE), axis.text.y = element_text(size = AXES_SIZE, face = "italic")) +
    # Now explicitly uses the title variable passed to the function
    labs(title = RESOLVED_TITLE$all, subtitle = subtitle, x = paste("Total", toupper(target_col)), y = RANK) +
    scale_fill_manual(values = pal_global$colours, drop = FALSE)

  # --- Plot 2: Utilisation propre du Helper stacked_bar avec facettes ! ---
  plot2_dt <- dt_merged[!is.na(name) & !is.na(date), .(Total_Abund = sum(get(target_col), na.rm = TRUE)), by = c("name", "date", RANK)]
  plot2_dt[, Taxa_Grouped := fifelse(get(RANK) %in% top_genera, get(RANK), "Others")]
  plot2_dt <- plot2_dt[, .(Total_Abund = sum(Total_Abund)), by = .(name, date, Taxa_Grouped)]

  plot2_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = taxa_order)]
  # Passed title and subtitle variables safely downstream
  p_stacked <- stacked_bar(plot2_dt, "date", "Total_Abund", "Taxa_Grouped", pal_global, RESOLVED_TITLE$each, subtitle, "Sampling Date", "Abundance", RANK, facet_col = "name")

  pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
  on.exit(if (names(dev.cur()) != "null device") dev.off())
  print(p_global); print(p_stacked)
}

# ==========================================================================
# Dispatch
# ==========================================================================
switch(MODE,
  "pathway_relative" = run_pathway_relative(),
  "relative_by_sample" = run_relative_by_sample(),
  "absolute_global" = run_absolute_global(),
  stop("Unknown mode: ", MODE)
)

message("✓ Execution completed. Output written to ", PDF, " and ", PARQUET)
