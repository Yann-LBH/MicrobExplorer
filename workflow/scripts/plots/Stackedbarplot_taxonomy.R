################################################################################
# Project : "MicrobExplorer"
# Script: "Heatmap"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(ggplot2)
library(arrow)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================
try(snakemake, silent = TRUE)
if (!exists("snakemake")) {
  snakemake <- list(
    input  = list(
      intersec = "Statistics/52.Final_Intersection_TaxaName",
      contigs  = "Statistics/3.RPKM_Contigs"
    ),
    output = list(
      pdf     = "Graphique/Barplot/Barplots.pdf",
      parquet = "Data/Parquet/barplot_taxonomy.parquet"
    ),
    params = list(
      target_rank = "Phylum",
      top_n       = 10L
    )
  )
}

setDTthreads(0L)

target_rank <- snakemake$params$target_rank
top_n       <- as.integer(snakemake$params$top_n)

TAX_RANKS <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

dir.create(dirname(snakemake$output$pdf),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake$output$parquet), recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 1. Chargement — une seule fonction, deux sources
# ==========================================================================
load_source <- function(path, source_label) {
  files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0L) {
    warning("Aucun fichier dans : ", path)
    return(NULL)
  }
  rbindlist(lapply(files, function(f) {
    dt        <- fread(f, sep = "\t", showProgress = FALSE)
    parts     <- strsplit(gsub("\\.tsv$", "", basename(f)), "_")[[1]]
    dt[, `:=`(
      Original_Taxon_Label = names(dt)[1],
      Source_Type          = source_label,
      Digesteur            = parts[3],
      Date_Raw             = parts[4],
      Date                 = as.Date(parts[4], format = "%y%m%d")
    )]
  }), use.names = TRUE, fill = TRUE)
}

sources <- list(
  "Graph Intersec" = snakemake$input$intersec,
  "Graph 3kb"      = snakemake$input$contigs
)

dt_all <- rbindlist(
  Filter(Negate(is.null), mapply(load_source, sources, names(sources), SIMPLIFY = FALSE)),
  use.names = TRUE, fill = TRUE
)

# ==========================================================================
# 2. Séparation taxonomie — vectorisée
# ==========================================================================
tax_split <- dt_all[, tstrsplit(Taxonomy, ";\\s*", fixed = FALSE,
                                names = TAX_RANKS, fill = "Unclassified")]
# Remplace vides par "Unclassified"
for (col in TAX_RANKS) {
  tax_split[get(col) == "" | get(col) == " ", (col) := "Unclassified"]
}
dt_taxo <- cbind(dt_all, tax_split)

# ==========================================================================
# 3. Parquet — format long annoté pour Shiny
# ==========================================================================
write_parquet(dt_taxo, snakemake$output$parquet)

# ==========================================================================
# 4. Agrégation
# ==========================================================================
setnames(dt_taxo, target_rank, "Taxon")

# Somme RPKM par échantillon x taxon
agg <- dt_taxo[, .(RPKM_Sum = sum(RPKM, na.rm = TRUE)),
               by = .(Source_Type, Digesteur, Date, Date_Raw, Taxon, Original_Taxon_Label)]

# Abondance relative par échantillon
agg[, Abundance_Pct := (RPKM_Sum / sum(RPKM_Sum)) * 100, by = .(Digesteur, Date)]

# Top N global par source
top_taxa <- agg[, .(Global_Sum = sum(RPKM_Sum)), by = .(Source_Type, Taxon)][
  , .SD[order(-Global_Sum)[seq_len(min(top_n, .N))]], by = Source_Type
][, .(Source_Type, Taxon)]

agg[top_taxa, Taxon_Final := Taxon, on = .(Source_Type, Taxon)]
agg[is.na(Taxon_Final), Taxon_Final := "Others"]

# Ré-agrégation après regroupement Others
final_dt <- agg[, .(Abundance_Pct = sum(Abundance_Pct)),
                by = .(Source_Type, Digesteur, Date, Date_Raw, Taxon_Final, Original_Taxon_Label)]

# Pré-split pour éviter les filtres dans la boucle
split_data <- split(final_dt, by = c("Source_Type", "Digesteur"), flatten = FALSE)

# ==========================================================================
# 5. Graphiques -> PDF compilé
# ==========================================================================
pdf(snakemake$output$pdf, width = 12, height = 8)

lapply(names(split_data), function(src) {
  lapply(names(split_data[[src]]), function(r) {
    
    plot_dt <- split_data[[src]][[r]]
    if (nrow(plot_dt) == 0L) return(invisible(NULL))
    
    # Ordre légende : abondance décroissante, Others en dernier
    taxon_order <- plot_dt[, .(total = sum(Abundance_Pct)), by = Taxon_Final][
      order(-total), Taxon_Final]
    taxon_order <- c(setdiff(taxon_order, "Others"),
                     intersect("Others", taxon_order))
    plot_dt[, Taxon_Final := factor(Taxon_Final, levels = taxon_order)]
    setorder(plot_dt, Date)
    
    p <- ggplot(plot_dt,
                aes(x = as.factor(Date), y = Abundance_Pct, fill = Taxon_Final)) +
      geom_col(position = "stack") +
      scale_fill_viridis_d(option = "turbo") +
      labs(
        title = paste("Abondance :", target_rank, "| Digesteur", r, "|", src),
        x     = "Date",
        y     = "Abondance Relative (%)",
        fill  = target_rank
      ) +
      theme_classic() +
      theme(
        plot.title    = element_text(size = 16, face = "bold"),
        axis.title    = element_text(size = 12, face = "bold"),
        axis.text     = element_text(size = 10),
        legend.title  = element_text(size = 12, face = "bold"),
        legend.text   = element_text(size = 10)
      )
    
    print(p)
  })
})

dev.off()

message("✓ PDF     : ", snakemake$output$pdf)
message("✓ Parquet : ", snakemake$output$parquet)