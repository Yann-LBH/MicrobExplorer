################################################################################
# Project : "MicrobExplorer"
# Script  : "Unified stacked barplot — pathway / taxonomy / organisms abundance"
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
library(data.table)
library(readxl)
library(ggplot2)
library(rlang)
library(arrow)
library(glue)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================

# Inputs
DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])[1]

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
TITLE_TEMPLATE    <- as.character(snakemake@params[["title"]])[1] %||% "{source} | Abundance of the top {top_n} {feature_type} in each sample"
SUBTITLE_TEMPLATE <- as.character(snakemake@params[["subtitle"]])[1] %||% "Mode : {mode} | Metric: {stand_col} | {rank}"
MODE              <- as.character(snakemake@params[["mode"]])[1]
TOP_N             <- as.integer(snakemake@params[["top_n"]])[1] %||% 10
STAND_COL         <- tolower(as.character(snakemake@params[["stand_col"]])[1])
RANK              <- as.character(snakemake@params[["rank"]])[1]

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# Standardization of Structure Type
FEATURE_TYPE <- if (SOURCE == "kegg") "pathways" else "taxon"

# TITRES ET SUBTITLES UNIFIÉS
RESOLVED_TITLE    <- glue(TITLE_TEMPLATE, source = toupper(SOURCE), top_n = TOP_N, feature_type = FEATURE_TYPE)
RESOLVED_SUBTITLE <- glue(SUBTITLE_TEMPLATE, mode = MODE, stand_col = toupper(STAND_COL), rank = RANK)

color_palette <- function(categories) {
  top_cats <- sort(setdiff(categories, "Others"))
  lvl_order <- c(top_cats, "Others")
  colours <- c(viridisLite::viridis(length(top_cats), option = tolower(PALETTE)), "#000000")
  names(colours) <- lvl_order
  list(colours = colours, levels = lvl_order)
}

# Stacked bar universel (Gère les facettes si demandées)
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
  
  # Si on passe un nom de couleur customisé ou généré
  if (is.list(colours)) colours <- colours$colours
  p <- p + scale_fill_manual(values = colours)
  
  # Ajout dynamique des facettes (Utilisé dans le mode Reads)
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(as.formula(paste0("~", facet_col)), scales = "free_x") +
      theme(strip.text = element_text(face = "bold", size = 12), panel.spacing = unit(1, "lines"))
  }
  return(p)
}

load_tsv_dir_dynamic <- function(paths, meta_dt, select_cols = NULL) {
  rbindlist(lapply(paths, function(f) {
    file_name <- basename(f)
    matched_sample <- meta_dt[sapply(sample_id, function(sid) grepl(sid, file_name)), sample_id]
    if (length(matched_sample) == 0 || is.na(matched_sample)) return(NULL)
    dt <- if (is.null(select_cols)) fread(f) else fread(f, select = intersect(select_cols, names(fread(f, nrows = 0))))
    if (nrow(dt) == 0) return(NULL)
    dt[, sample_id := matched_sample]
  }), use.names = TRUE, fill = TRUE)
}

load_metadata <- function(path) {
  dt <- if (tolower(tools::file_ext(path)) %in% c("xlsx", "xls")) as.data.table(readxl::read_excel(path)) else fread(path)
  dt[, c("sample_id", "name", "date") := .(as.character(sample_id), as.character(name), as.character(date))]
}

# ==========================================================================
# RUN MODES
# ==========================================================================

run_pathway_relative <- function(title = RESOLVED_TITLE, subtitle = RESOLVED_SUBTITLE) {
  target_col <- STAND_COL
  if (is.na(target_col) || target_col == "") stop("L'argument 'stand_col' est manquant ou égal à NA.")

  meta <- load_metadata(METADATA)
  all_data <- load_tsv_dir_dynamic(DATA, meta, select_cols = c(RANK, target_col))
  write_parquet(all_data[meta, on = "sample_id", nomatch = 0L], PARQUET)

  agg <- all_data[meta, on = "sample_id", nomatch = 0L][, .(total = sum(get(target_col), na.rm = TRUE)), by = .(name, date, get(RANK))]
  setnames(agg, "get", RANK)
  agg[, category := fifelse(frank(-total, ties.method = "random") <= TOP_N, as.character(get(RANK)), "Others"), by = .(name, date)]
  
  final_dt <- agg[, .(sum_val = sum(total)), by = .(name, date, category)][, pct := (sum_val / sum(sum_val)) * 100, by = .(name, date)]
  pal <- color_palette(unique(final_dt$category))
  final_dt[, category := factor(category, levels = pal$levels)]

  pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
  lapply(split(final_dt, by = "name"), function(df_cond) {
    # Using dynamic title matching the exact sample name
    print(stacked_bar(df_cond, "date", "pct", "category", pal, title, subtitle, "Date", "Relative Abundance (%)", "Pathways"))
  })
  dev.off()
}

run_relative_by_sample <- function(title = RESOLVED_TITLE, subtitle = RESOLVED_SUBTITLE) {
  target_col <- STAND_COL
  if (is.na(target_col) || target_col == "") stop("L'argument 'stand_col' est manquant ou égal à NA.")

  meta <- load_metadata(METADATA)
  dt_taxo <- load_tsv_dir_dynamic(DATA, meta)[meta, on = "sample_id", nomatch = 0L]
  write_parquet(dt_taxo, PARQUET)

  setnames(dt_taxo, RANK, "Taxon")
  agg <- dt_taxo[, .(Abund_Sum = sum(get(target_col), na.rm = TRUE)), by = .(name, date, Taxon)][, Abund_Pct := (Abund_Sum / sum(Abund_Sum)) * 100, by = .(name, date)]
  top_taxa <- agg[, .(G = sum(Abund_Sum)), by = Taxon][order(-G)[seq_len(min(TOP_N, .N))], Taxon]
  agg[, Taxon_Final := fifelse(Taxon %in% top_taxa, Taxon, "Others")]
  final_dt <- agg[, .(Abund_Pct = sum(Abund_Pct)), by = .(name, date, Taxon_Final)]

  pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
  # Harmonized: loop using split() instead of unique() + manual filtering
  lapply(split(final_dt, by = "name"), function(df_cond) {
    pal <- color_palette(unique(df_cond$Taxon_Final))
    df_cond[, Taxon_Final := factor(Taxon_Final, levels = pal$levels)]
    # Standardized title structure to match the sample view
    print(stacked_bar(df_cond, "date", "Abund_Pct", "Taxon_Final", pal, title, subtitle, "Date", "Relative Abundance (%)", RANK))
  })
  dev.off()
}

run_absolute_global <- function(title = RESOLVED_TITLE, subtitle = RESOLVED_SUBTITLE) {
  target_col <- STAND_COL
  if (is.na(target_col) || target_col == "") stop("L'argument 'stand_col' est manquant ou égal à NA.")

  meta <- load_metadata(METADATA)
  dt <- load_tsv_dir_dynamic(DATA, meta)[meta, on = "sample_id", nomatch = 0L]
  if (!RANK %in% names(dt)) stop(sprintf("The taxonomy column [%s] is missing.", RANK))
  dt[get(RANK) == "" | is.na(get(RANK)), (RANK) := "Unclassified"]
  write_parquet(dt, PARQUET)

  taxa_config_global <- dt[, .(Global_Abund = sum(get(target_col), na.rm = TRUE)), by = c(RANK)][order(-Global_Abund)]
  top_genera <- taxa_config_global[seq_len(min(TOP_N, .N))][[RANK]]

  # --- Plot 1: Global horizontal bar ---
  plot1_dt <- taxa_config_global[, .(Taxa_Grouped = fifelse(get(RANK) %in% top_genera, get(RANK), "Others"), Global_Abund)][, .(Total_Abund = sum(Global_Abund)), by = Taxa_Grouped][order(Total_Abund)]
  plot1_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = Taxa_Grouped)]
  
  p_global <- ggplot(plot1_dt, aes(x = Total_Abund, y = Taxa_Grouped, fill = Taxa_Grouped)) +
    geom_col(show.legend = FALSE) + get(THEME)() + theme(plot.title = element_text(face = "bold", size = TITLE_SIZE), axis.text.y = element_text(size = AXES_SIZE, face = "italic")) +
    # Now explicitly uses the title variable passed to the function
    labs(title = title, subtitle = subtitle, x = paste("Total", toupper(target_col)), y = RANK) +
    scale_fill_manual(values = color_palette(unique(plot1_dt$Taxa_Grouped))$colours)

  # --- Plot 2: Utilisation propre du Helper stacked_bar avec facettes ! ---
  plot2_dt <- dt[!is.na(name) & !is.na(date), .(Total_Abund = sum(get(target_col), na.rm = TRUE)), by = c("name", "date", RANK)]
  plot2_dt[, Taxa_Grouped := fifelse(get(RANK) %in% top_genera, get(RANK), "Others")]
  plot2_dt <- plot2_dt[, .(Total_Abund = sum(Total_Abund)), by = .(name, date, Taxa_Grouped)]
  pal_stacked <- color_palette(unique(plot2_dt$Taxa_Grouped))
  plot2_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = pal_stacked$levels)]

  # Passed title and subtitle variables safely downstream
  p_stacked <- stacked_bar(plot2_dt, "date", "Total_Abund", "Taxa_Grouped", pal_stacked, title, subtitle, "Sampling Date", "Abundance", RANK, facet_col = "name")

  pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
  print(p_global); print(p_stacked)
  dev.off()
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
