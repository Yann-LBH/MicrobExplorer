################################################################################
# Project : "MicrobExplorer"
# Script: "Utilitaire : Barplot_mapping_percent"
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
    input  = list(txt_dir = "Statistics/Coassembly"),
    output = list(
      pdf     = "Graphique/Barplot/Statistics/Stats_mapping.pdf",
      parquet = "Data/Parquet/mapping_rates.parquet"
    ),
    params = list(
      pattern = "Pourcentage de mapping"
    )
  )
}

setDTthreads(0L)

dir.create(dirname(snakemake$output$pdf),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake$output$parquet), recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 1. Chargement
# ==========================================================================
files <- list.files(snakemake$input$txt_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(files) == 0L) stop("Aucun fichier trouvé dans : ", snakemake$input$txt_dir)

dt <- rbindlist(lapply(files, function(f) {
  lines      <- readLines(f, warn = FALSE)
  sample     <- gsub("\\.txt$", "", basename(f))
  sample     <- gsub("counted_mapped_", "", sample)
  sample     <- gsub("_([0-9]{2})([0-9]{2})([0-9]{2})", " 20\\1-\\2-\\3", sample)
  cible      <- grep(snakemake$params$pattern, lines, value = TRUE)
  valeur     <- as.numeric(gsub(",", ".", gsub("[^0-9.,]", "", cible)))
  data.table(Echantillon = sample, Mapping_Rate = valeur)
}))

# Ordre décroissant pour le plot
dt[, Echantillon := factor(Echantillon, levels = rev(Echantillon))]

# ==========================================================================
# 2. Parquet
# ==========================================================================
write_parquet(dt, snakemake$output$parquet)

# ==========================================================================
# 3. Graphique -> PDF
# ==========================================================================
p <- ggplot(dt, aes(x = Echantillon, y = Mapping_Rate)) +
  geom_col(fill = "#008c8e") +
  coord_flip() +
  labs(
    title = "Taux de Mapping par Échantillon",
    x     = NULL,
    y     = "Pourcentage (%)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

pdf(snakemake$output$pdf, width = 10, height = 6)
print(p)
dev.off()

message("✓ PDF     : ", snakemake$output$pdf)
message("✓ Parquet : ", snakemake$output$parquet)