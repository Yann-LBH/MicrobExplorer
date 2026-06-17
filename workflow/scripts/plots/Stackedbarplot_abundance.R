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

library(data.table)
library(ggplot2)
library(arrow)

# ==========================================================================
# Snakemake configuration
# ==========================================================================

# Inputs
DATA <- as.character(snakemake@input$data)
METADATA <- as.character(snakemake@input$metadata)

# Outputs
PDF <- as.character(snakemake@output$pdf)
PARQUET <- as.character(snakemake@output$parquet)

# Parameters
MODE <- as.character(snakemake@params$mode)
TOP_N <- as.integer(snakemake@params$top_n)
VALUE_COL <- as.character(snakemake@params$value_col)
TARGET_RANK <- as.character(snakemake@params$target_rank) # taxonomy mode
TAXON_RANK <- as.integer(snakemake@params$taxon_rank) # organisms mode

TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

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

# Generic Top-N collapse: returns a data.table with a `category` column
top_n_collapse <- function(dt, value_col, group_col, by_cols, top_n) {
  # Aggregate totals per group to rank
  agg <- dt[, .(total = sum(get(value_col), na.rm = TRUE)), by = c(by_cols, group_col)]
  agg[, category := {
    r <- frank(-total, ties.method = "random")
    fifelse(r <= top_n, get(group_col), "Others")
  }, by = by_cols]
  agg[, .(sum_val = sum(total)), by = c(by_cols, "category")]
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

load_tsv_dir <- function(paths, select_cols = NULL) {
  rbindlist(
    lapply(paths, function(f) {
      dt <- if (is.null(select_cols)) {
        fread(f, showProgress = FALSE)
      } else {
        fread(f, select = select_cols, showProgress = FALSE)
      }
      dt[, sample_id := gsub("^[^_]+_|\\.tsv$", "", basename(f))]
      dt
    }),
    use.names = TRUE, fill = TRUE
  )
}

load_metadata <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    as.data.table(readxl::read_excel(path))
  } else {
    fread(path)
  }
}

# ==========================================================================
# MODE: pathway
# ==========================================================================
run_pathway <- function() {
  files <- unlist(snakemake$input$data)
  meta <- load_metadata(snakemake$input$metadata)

  all_data <- load_tsv_dir(files, select_cols = c("level_3", VALUE_COL))

  setkey(all_data, sample_id)
  setkey(meta, sample_id)
  all_data <- all_data[meta, nomatch = 0L]

  write_parquet(all_data, snakemake$output$parquet)

  # Aggregate per condition × date × pathway
  agg <- all_data[, .(total = sum(get(VALUE_COL), na.rm = TRUE)),
    by = .(name, date, level_3)
  ]

  agg[, category := {
    r <- frank(-total, ties.method = "random")
    fifelse(r <= TOP_N, level_3, "Others")
  }, by = .(name, date)]

  final_dt <- agg[, .(sum_val = sum(total)), by = .(name, date, category)]
  final_dt[, pct := (sum_val / sum(sum_val)) * 100, by = .(name, date)]

  pal <- make_palette(unique(final_dt$category))
  final_dt[, category := factor(category, levels = pal$levels)]

  pdf(snakemake$output$pdf, width = 12, height = 8)
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
# MODE: taxonomy
# ==========================================================================
run_taxonomy <- function() {
  load_source <- function(path, source_label) {
    files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)
    if (!length(files)) {
      warning("No files in: ", path)
      return(NULL)
    }
    rbindlist(lapply(files, function(f) {
      dt <- fread(f, sep = "\t", showProgress = FALSE)
      parts <- strsplit(gsub("\\.tsv$", "", basename(f)), "_")[[1]]
      dt[, `:=`(
        Original_Taxon_Label = names(dt)[1],
        Source_Type          = source_label,
        Digesteur            = parts[3],
        Date                 = as.Date(parts[4], format = "%y%m%d")
      )]
    }), use.names = TRUE, fill = TRUE)
  }

  sources <- list(
    "Graph Intersec" = snakemake$input$intersec,
    "Graph 3kb"      = snakemake$input$contigs
  )

  dt_all <- rbindlist(
    Filter(
      Negate(is.null),
      mapply(load_source, sources, names(sources), SIMPLIFY = FALSE)
    ),
    use.names = TRUE, fill = TRUE
  )

  # Vectorised taxonomy split
  tax_split <- dt_all[, tstrsplit(Taxonomy, ";\\s*",
    fixed = FALSE,
    names = TAX_RANKS, fill = "Unclassified"
  )]
  for (col in TAX_RANKS) {
    tax_split[get(col) == "" | get(col) == " ", (col) := "Unclassified"]
  }

  dt_taxo <- cbind(dt_all, tax_split)
  write_parquet(dt_taxo, snakemake$output$parquet)

  setnames(dt_taxo, TARGET_RANK, "Taxon")

  agg <- dt_taxo[, .(RPKM_Sum = sum(RPKM, na.rm = TRUE)),
    by = .(Source_Type, Digesteur, Date, Taxon, Original_Taxon_Label)
  ]
  agg[, Abund_Pct := (RPKM_Sum / sum(RPKM_Sum)) * 100, by = .(Digesteur, Date)]

  # Top N per source
  top_taxa <- agg[, .(G = sum(RPKM_Sum)), by = .(Source_Type, Taxon)][
    , .SD[order(-G)[seq_len(min(TOP_N, .N))]],
    by = Source_Type
  ][, .(Source_Type, Taxon)]

  agg[top_taxa, Taxon_Final := Taxon, on = .(Source_Type, Taxon)]
  agg[is.na(Taxon_Final), Taxon_Final := "Others"]

  final_dt <- agg[, .(Abund_Pct = sum(Abund_Pct)),
    by = .(Source_Type, Digesteur, Date, Taxon_Final, Original_Taxon_Label)
  ]

  pdf(snakemake$output$pdf, width = 12, height = 8)
  lapply(unique(final_dt$Source_Type), function(src) {
    lapply(unique(final_dt[Source_Type == src, Digesteur]), function(r) {
      plot_dt <- final_dt[Source_Type == src & Digesteur == r]
      if (!nrow(plot_dt)) {
        return(invisible(NULL))
      }

      taxon_order <- c(
        setdiff(plot_dt[, .(t = sum(Abund_Pct)), by = Taxon_Final][order(-t), Taxon_Final], "Others"),
        "Others"
      )
      plot_dt[, Taxon_Final := factor(Taxon_Final, levels = taxon_order)]
      setorder(plot_dt, Date)

      pal <- setNames(
        c(viridisLite::viridis(length(taxon_order) - 1L, option = "turbo"), "#000000"),
        taxon_order
      )

      print(stacked_bar(plot_dt, "Date", "Abund_Pct", "Taxon_Final",
        pal,
        title = paste("Abundance:", TARGET_RANK, "| Digesteur", r, "|", src),
        subtitle = NULL,
        x_lab = "Date", y_lab = "Relative Abundance (%)",
        fill_lab = TARGET_RANK
      ))
    })
  })
  dev.off()
}

# ==========================================================================
# MODE: organisms
# ==========================================================================
run_organisms <- function() {
  files <- unlist(snakemake$input$data)
  meta <- load_metadata(snakemake$input$metadata)

  dt <- load_tsv_dir(files)

  # Vectorised genus extraction
  dt[, Genus := {
    g <- trimws(sapply(strsplit(Taxonomy, ";"), `[`, TAXON_RANK))
    fifelse(is.na(g) | g %in% c("NA", ""), "Unclassified Genus", g)
  }]

  setkey(dt, sample_id)
  setkey(meta, sample_id)
  dt <- meta[dt, nomatch = 0L]

  write_parquet(dt, snakemake$output$parquet)

  # Global top genera
  genus_global <- dt[, .(Global_RPKM = sum(RPKM, na.rm = TRUE)), by = Genus]
  setorder(genus_global, -Global_RPKM)
  top_genera <- genus_global[seq_len(min(TOP_N, .N)), Genus]

  # --- Plot 1: global horizontal bar ---
  plot1_dt <- genus_global[, .(
    Genus_Grouped = fifelse(Genus %in% top_genera, Genus, "Others"),
    Global_RPKM
  )][, .(Total_RPKM = sum(Global_RPKM)), by = Genus_Grouped]
  setorder(plot1_dt, Total_RPKM)
  plot1_dt[, Genus_Grouped := factor(Genus_Grouped, levels = Genus_Grouped)]

  p_global <- ggplot(
    plot1_dt,
    aes(x = Total_RPKM, y = Genus_Grouped, fill = Genus_Grouped)
  ) +
    geom_col(show.legend = FALSE) +
    scale_fill_viridis_d(option = "turbo") +
    labs(
      title = sprintf("Global Abundance — Top %d Genera", TOP_N),
      subtitle = "Aggregated data from all digesters",
      x = "Total RPKM (Summed)", y = "Genus"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, face = "italic"),
      plot.title = element_text(face = "bold", size = 14)
    )

  # --- Plot 2: stacked per digester × date ---
  plot2_dt <- dt[!is.na(name) & !is.na(date),
    .(Total_RPKM = sum(RPKM, na.rm = TRUE)),
    by = .(name, date, Genus)
  ]
  plot2_dt[, Genus_Grouped := fifelse(Genus %in% top_genera, Genus, "Others")]
  plot2_dt <- plot2_dt[, .(Total_RPKM = sum(Total_RPKM)), by = .(name, date, Genus_Grouped)]
  setorder(plot2_dt, Total_RPKM)
  plot2_dt[, Genus_Grouped := factor(Genus_Grouped, levels = unique(Genus_Grouped))]

  p_stacked <- ggplot(
    plot2_dt,
    aes(x = as.factor(date), y = Total_RPKM, fill = Genus_Grouped)
  ) +
    geom_col(colour = "white", linewidth = 0.1) +
    scale_fill_viridis_d(option = "turbo") +
    facet_wrap(~name, scales = "free_x") +
    labs(
      title = "Composition of Bacterial Genera by Digester",
      subtitle = "Vertical stacked bars | RPKM values",
      x = "Sampling Date", y = "Total RPKM", fill = "Genera"
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

  pdf(snakemake$output$pdf, width = 14, height = 8)
  print(p_global)
  print(p_stacked)
  dev.off()
}

# ==========================================================================
# Dispatch
# ==========================================================================
switch(MODE,
  pathway   = run_pathway(),
  taxonomy  = run_taxonomy(),
  organisms = run_organisms(),
  stop("Unknown mode: ", MODE, ". Use 'pathway', 'taxonomy', or 'organisms'.")
)

message("Done. PDF: ", snakemake$output$pdf)
message("Done. Parquet: ", snakemake$output$parquet)
