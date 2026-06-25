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

# All script comments are provided in English as requested.
library(data.table)
library(readxl)
library(ggplot2)
library(arrow)

# Inputs
DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])

# Outputs
PDF     <- as.character(c(snakemake@output[["pdf"]])[1])
PARQUET <- as.character(c(snakemake@output[["parquet"]])[1])

# Parameters
MODE          <- as.character(c(snakemake@params[["mode"]])[1])
TOP_N         <- as.integer(c(snakemake@params[["top_n"]])[1])
VALUE_COL     <- as.character(c(snakemake@params[["value_col"]])[1])
TARGET_RANK   <- tolower(as.character(c(snakemake@params[["target_rank"]])[1])) %||% "genus"
TAXON_RANK    <- tolower(as.character(c(snakemake@params[["taxon_rank"]])[1])) %||% "genus"
PATHWAY_LEVEL <- as.character(c(snakemake@params[["pathway_level"]])[1])

# ✅ FIXED: Define TAX_RANKS globally so it is accessible in all modes
TAX_RANKS     <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# ==========================================================================
# Helpers
# ==========================================================================

# Build turbo palette: top categories in colour, "Others" in black
make_palette <- function(categories) {
  top_cats <- sort(setdiff(categories, "Others"))
  lvl_order <- c(top_cats, "Others")
  colours <- c(viridisLite::viridis(length(top_cats), option = "turbo"), "#000000")
  names(colours) <- lvl_order
  list(colours = colours, levels = lvl_order)
}

# Stacked bar ggplot
stacked_bar <- function(dt, x_col, y_col, fill_col, colours, title, subtitle,
                        x_lab, y_lab, fill_lab) {
  ggplot(dt, aes(x = as.factor(get(x_col)), y = get(y_col), fill = get(fill_col))) +
    geom_bar(
      stat = "identity",
      position = position_stack(reverse = TRUE),
      colour = "white",
      linewidth = 0.05
    ) +
    scale_fill_manual(values = colours) +
    labs(
      title = title, subtitle = subtitle,
      x = x_lab, y = y_lab, fill = fill_lab
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 7),
      panel.grid.major.x = element_blank()
    )
}

# ==========================================================================
# Load helpers — shared across modes
# ==========================================================================

load_tsv_dir_dynamic <- function(paths, meta_dt, select_cols = NULL) {
  rbindlist(
    lapply(paths, function(f) {
      file_name <- basename(f)
      matched_sample <- meta_dt[sapply(sample_id, function(sid) grepl(sid, file_name)), sample_id]
      
      if (length(matched_sample) == 0 || is.na(matched_sample)) return(NULL)
      
      dt <- if (is.null(select_cols)) {
        fread(f, showProgress = FALSE)
      } else {
        actual_cols <- intersect(select_cols, names(fread(f, nrows = 0)))
        fread(f, select = actual_cols, showProgress = FALSE)
      }
      
      if (nrow(dt) == 0) return(NULL)
      
      dt[, sample_id := matched_sample]
      dt
    }),
    use.names = TRUE, fill = TRUE
  )
}

load_metadata <- function(path) {
  ext <- tolower(tools::file_ext(path))
  dt  <- if (ext %in% c("xlsx", "xls")) as.data.table(readxl::read_excel(path)) else fread(path)
  dt[, sample_id := as.character(sample_id)]
  dt[, name      := as.character(name)]
  dt[, date      := as.character(date)]
  dt
}

# ==========================================================================
# MODE: Kegg
# ==========================================================================
run_kegg <- function() {
  meta <- load_metadata(METADATA)
  all_data <- load_tsv_dir_dynamic(DATA, meta, select_cols = c(PATHWAY_LEVEL, VALUE_COL))

  setkey(all_data, sample_id)
  setkey(meta, sample_id)
  all_data <- all_data[meta, nomatch = 0L]

  write_parquet(all_data, PARQUET)

  agg <- all_data[, .(total = sum(get(VALUE_COL), na.rm = TRUE)),
    by = .(name, date, get(PATHWAY_LEVEL))
  ]
  setnames(agg, "get", PATHWAY_LEVEL)

  agg[, category := {
    r <- frank(-total, ties.method = "random")
    fifelse(r <= TOP_N, as.character(get(PATHWAY_LEVEL)), "Others")
  }, by = .(name, date)]

  final_dt <- agg[, .(sum_val = sum(total)), by = .(name, date, category)]
  final_dt[, pct := (sum_val / sum(sum_val)) * 100, by = .(name, date)]

  pal <- make_palette(unique(final_dt$category))
  final_dt[, category := factor(category, levels = pal$levels)]

  pdf(PDF, width = 12, height = 8)
  lapply(split(final_dt, by = "name", keep.by = TRUE), function(df_cond) {
    print(stacked_bar(df_cond, "date", "pct", "category",
      pal$colours,
      title = paste("Condition:", df_cond$name[1L]),
      subtitle = paste("Top", TOP_N, "abundance pathways"),
      x_lab = "Date", y_lab = "Relative Abundance (%)",
      fill_lab = "Pathways"
    ))
  })
  dev.off()
}

# ==========================================================================
# MODE: Contigs
# ==========================================================================
run_contigs <- function() {
  meta   <- load_metadata(METADATA)
  dt_all <- load_tsv_dir_dynamic(DATA, meta)

  tax_split <- dt_all[, tstrsplit(Taxonomy, ";\\s*",
    fixed = FALSE,
    names = TAX_RANKS, fill = "Unclassified"
  )]
  
  for (col in TAX_RANKS) {
    tax_split[get(col) == "" | get(col) == " " | is.na(get(col)), (col) := "Unclassified"]
  }

  dt_taxo <- cbind(dt_all, tax_split)
  
  setkey(dt_taxo, sample_id)
  setkey(meta, sample_id)
  dt_taxo <- meta[dt_taxo, nomatch = 0L]
  
  write_parquet(dt_taxo, PARQUET)

  setnames(dt_taxo, TARGET_RANK, "Taxon")

  agg <- dt_taxo[, .(RPKM_Sum = sum(get(VALUE_COL), na.rm = TRUE)),
    by = .(name, date, Taxon)
  ]
  agg[, Abund_Pct := (RPKM_Sum / sum(RPKM_Sum)) * 100, by = .(name, date)]

  top_taxa <- agg[, .(G = sum(RPKM_Sum)), by = Taxon][
    order(-G)[seq_len(min(TOP_N, .N))], Taxon
  ]

  agg[, Taxon_Final := fifelse(Taxon %in% top_taxa, Taxon, "Others")]

  final_dt <- agg[, .(Abund_Pct = sum(Abund_Pct)),
    by = .(name, date, Taxon_Final)
  ]

  pdf(PDF, width = 12, height = 8)
  
  lapply(unique(final_dt$name), function(r) {
    plot_dt <- final_dt[name == r]
    if (!nrow(plot_dt)) {
      return(invisible(NULL))
    }

    taxon_order <- c(
      setdiff(plot_dt[, .(t = sum(Abund_Pct)), by = Taxon_Final][order(-t), Taxon_Final], "Others"),
      "Others"
    )
    plot_dt[, Taxon_Final := factor(Taxon_Final, levels = taxon_order)]
    setorder(plot_dt, date)

    pal <- setNames(
      c(viridisLite::viridis(length(taxon_order) - 1L, option = "turbo"), "#000000"),
      taxon_order
    )

    print(stacked_bar(plot_dt, "date", "Abund_Pct", "Taxon_Final",
      pal,
      title = paste("Abundance:", TARGET_RANK, "| Digesteur", r),
      subtitle = NULL,
      x_lab = "Date", y_lab = "Relative Abundance (%)",
      fill_lab = TARGET_RANK
    ))
  })
  
  dev.off()
}

# ==========================================================================
# MODE: Reads
# ==========================================================================
run_reads <- function() {
  meta <- load_metadata(METADATA)
  dt   <- load_tsv_dir_dynamic(DATA, meta)

  if (!TAXON_RANK %in% names(dt)) {
    stop(sprintf("The taxonomy column [%s] is missing from the input file.", TAXON_RANK))
  }

  dt[get(TAXON_RANK) == "" | get(TAXON_RANK) == " " | is.na(get(TAXON_RANK)), (TAXON_RANK) := "Unclassified"]

  setkey(dt, sample_id)
  setkey(meta, sample_id)
  dt <- meta[dt, nomatch = 0L]

  write_parquet(dt, PARQUET)

  # ✅ FIXED: Standardized dynamic column extraction with data.table evaluating via character vector
  taxa_config_global <- dt[, .(Global_RPKM = sum(get(VALUE_COL), na.rm = TRUE)), by = c(TAXON_RANK)]
  setorder(taxa_config_global, -Global_RPKM)
  
  # Use standard vector indexing to prevent length/evaluation issues
  top_genera <- taxa_config_global[seq_len(min(TOP_N, .N)), [[TAXON_RANK]]]

  # --- Plot 1: global horizontal bar ---
  plot1_dt <- taxa_config_global[, .(
    Taxa_Grouped = fifelse(get(TAXON_RANK) %in% top_genera, get(TAXON_RANK), "Others"),
    Global_RPKM
  )][, .(Total_RPKM = sum(Global_RPKM)), by = Taxa_Grouped]
  setorder(plot1_dt, Total_RPKM)
  plot1_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = Taxa_Grouped)]

  p_global <- ggplot(
    plot1_dt,
    aes(x = Total_RPKM, y = Taxa_Grouped, fill = Taxa_Grouped)
  ) +
    geom_col(show.legend = FALSE) +
    scale_fill_viridis_d(option = "turbo") +
    labs(
      title = sprintf("Global Abundance — Top %d %ss", TOP_N, TAXON_RANK),
      subtitle = "Aggregated data from all matched samples",
      x = paste("Total", VALUE_COL, "(Summed)"), y = TAXON_RANK
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, face = "italic"),
      plot.title = element_text(face = "bold", size = 14)
    )

  # --- Plot 2: stacked per sample name × date ---
  # ✅ FIXED: Force evaluation context using character vector syntax
  plot2_dt <- dt[!is.na(name) & !is.na(date),
    .(Total_RPKM = sum(get(VALUE_COL), na.rm = TRUE)),
    by = c("name", "date", TAXON_RANK)
  ]
  
  plot2_dt[, Taxa_Grouped := fifelse(get(TAXON_RANK) %in% top_genera, get(TAXON_RANK), "Others")]
  plot2_dt <- plot2_dt[, .(Total_RPKM = sum(Total_RPKM)), by = .(name, date, Taxa_Grouped)]
  setorder(plot2_dt, Total_RPKM)
  plot2_dt[, Taxa_Grouped := factor(Taxa_Grouped, levels = unique(Taxa_Grouped))]

  p_stacked <- ggplot(
    plot2_dt,
    aes(x = as.factor(date), y = Total_RPKM, fill = Taxa_Grouped)
  ) +
    geom_col(colour = "white", linewidth = 0.1) +
    scale_fill_viridis_d(option = "turbo") +
    facet_wrap(~name, scales = "free_x") +
    labs(
      title = sprintf("Composition of Bacterial %ss", TAXON_RANK),
      subtitle = paste("Vertical stacked bars |", VALUE_COL, "values"),
      x = "Sampling Date", y = VALUE_COL, fill = TAXON_RANK
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 9, face = "italic"),
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      panel.spacing = unit(1, "lines"),
      plot.title = element_text(face = "bold", size = 16)
    )

  pdf(PDF, width = 14, height = 8)
  print(p_global)
  print(p_stacked)
  dev.off()
}

# ==========================================================================
# Dispatch
# ==========================================================================
switch(MODE,
  pathway   = run_kegg(),
  taxonomy  = run_contigs(),
  organisms = run_reads(),
  stop("Unknown mode: ", MODE, ". Use 'kegg', 'contigs', or 'reads'.")
)

message("✓ Execution completed. Output written to ", PDF, " and ", PARQUET)