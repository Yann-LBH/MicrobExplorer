################################################################################
# Project : "MicrobExplorer"
# Script: "Stackedbarplot organisms abundances"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(ggplot2)
library(arrow)
library(readxl)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================
try(snakemake, silent = TRUE)
if (!exists("snakemake")) {
  snakemake <- list(
    input  = list(
      tsv_dir  = "Statistics/6.Final_Results",
      metadata = "metadata/metadata.xlsx"
    ),
    output = list(
      pdf     = "Graphique/Barplot/Statistics/Barplots.pdf",
      parquet = "Data/Parquet/genus_abundance.parquet"
    ),
    params = list(
      top_n      = 19L,
      taxon_rank = 6L   # index du genre dans la chaine Taxonomy (";"-séparé)
    )
  )
}

setDTthreads(0L)

top_n      <- as.integer(snakemake$params$top_n)
taxon_rank <- as.integer(snakemake$params$taxon_rank)

dir.create(dirname(snakemake$output$pdf),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake$output$parquet), recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 1. Chargement — une seule passe pour les deux graphiques
# ==========================================================================
files <- list.files(snakemake$input$tsv_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0L) stop("Aucun fichier TSV trouvé dans : ", snakemake$input$tsv_dir)

meta <- as.data.table(read_excel(snakemake$input$metadata))

dt <- rbindlist(
  lapply(files, function(f) {
    d <- fread(f, sep = "\t", showProgress = FALSE)
    d[, sample_id := gsub("Final_Annotated_|\\.tsv$", "", basename(f))]
    d
  }),
  use.names = TRUE, fill = TRUE
)

# Extraction du genre vectorisée — pas de map_chr ligne par ligne
dt[, Genus := {
  g <- trimws(sapply(strsplit(Taxonomy, ";"), `[`, taxon_rank))
  fifelse(is.na(g) | g == "NA" | g == "", "Unclassified Genus", g)
}]

# Jointure avec metadata
setkey(dt,   sample_id)
setkey(meta, sample_id)
dt <- meta[dt, nomatch = 0L]

# ==========================================================================
# 2. Parquet — format long complet pour Shiny
# ==========================================================================
write_parquet(dt, snakemake$output$parquet)

# ==========================================================================
# 3. Top N global (partagé par les deux graphiques)
# ==========================================================================
genus_global <- dt[, .(Global_RPKM = sum(RPKM, na.rm = TRUE)), by = Genus]
setorder(genus_global, -Global_RPKM)
top_genera <- genus_global[seq_len(min(top_n, .N)), Genus]

# ==========================================================================
# 4. Données agrégées pour les plots
# ==========================================================================

# -- Plot 1 : Global bar (tous échantillons confondus)
plot1_dt <- genus_global[, .(
  Genus_Grouped = fifelse(Genus %in% top_genera, Genus, "Others"),
  Global_RPKM
)][, .(Total_RPKM = sum(Global_RPKM)), by = Genus_Grouped]
setorder(plot1_dt, Total_RPKM)
plot1_dt[, Genus_Grouped := factor(Genus_Grouped, levels = Genus_Grouped)]

# -- Plot 2 : Stacked par digesteur x date
plot2_dt <- dt[!is.na(name) & !is.na(date),
               .(Total_RPKM = sum(RPKM, na.rm = TRUE)),
               by = .(name, date, Genus)]
plot2_dt[, Genus_Grouped := fifelse(Genus %in% top_genera, Genus, "Others")]
plot2_dt <- plot2_dt[, .(Total_RPKM = sum(Total_RPKM)), by = .(name, date, Genus_Grouped)]
setorder(plot2_dt, Total_RPKM)
plot2_dt[, Genus_Grouped := factor(Genus_Grouped, levels = unique(Genus_Grouped))]

# ==========================================================================
# 5. Graphiques -> PDF compilé
# ==========================================================================
p_global <- ggplot(plot1_dt,
                   aes(x = Total_RPKM, y = Genus_Grouped, fill = Genus_Grouped)) +
  geom_col(show.legend = FALSE) +
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title    = sprintf("Global Abundance — Top %d Genera", top_n),
    subtitle = "Aggregated data from all digesters",
    x        = "Total RPKM (Summed)",
    y        = "Genus"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "italic"),
    plot.title  = element_text(face = "bold", size = 14)
  )

p_stacked <- ggplot(plot2_dt,
                    aes(x = as.factor(date), y = Total_RPKM, fill = Genus_Grouped)) +
  geom_col(color = "white", linewidth = 0.1) +
  scale_fill_viridis_d(option = "turbo") +
  facet_wrap(~ name, scales = "free_x") +
  labs(
    title    = "Composition of Bacterial Genera by Digester",
    subtitle = "Vertical stacked bars | RPKM values",
    x        = "Sampling Date",
    y        = "Total RPKM",
    fill     = "Genera"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 10),
    legend.position  = "right",
    legend.text      = element_text(size = 9, face = "italic"),
    legend.title     = element_text(face = "bold"),
    strip.text       = element_text(face = "bold", size = 12),
    panel.spacing    = unit(1, "lines"),
    plot.title       = element_text(face = "bold", size = 16)
  )

pdf(snakemake$output$pdf, width = 14, height = 8)
print(p_global)
print(p_stacked)
dev.off()

message("✓ PDF     : ", snakemake$output$pdf)
message("✓ Parquet : ", snakemake$output$parquet)