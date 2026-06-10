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
library(viridis)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================

#Inputs
DATA        <- snakemake[["input"]][["data"]]

#Ouputs
PDF <- snakemake[["output"]][["pdf"]]
PARQUET <- snakemake[["output"]][["parquet"]]

#Parameters
ACTIVE_MODULES <- snakemake[["params"]][["active_modules"]]

# ==========================================================================
# 1. Chargement
# ==========================================================================
dt_raw <- fread(DATA, sep = "|", strip.white = TRUE)

# ==========================================================================
# 2. Calcul des pertes — en place
# ==========================================================================
dt_raw[, `:=`(
  P_Counted        = Brut - Counted,
  P_Trimmed        = Counted - Trimmed,
  P_RPKM           = Counted - RPKM,
  P_RPKM_Filtered  = RPKM - RPKM_Filtered,
  Intersection     = RPKM_Filtered - Intersection
)]

dt_wide <- dt_raw[, .SD, .SDcols = c("Sample", ordre_final)]

# ==========================================================================
# 3. Parquet — format wide + long pour Shiny
# ==========================================================================
write_parquet(dt_wide, PARQUET)

# Format long pour le plot
dt_plot <- melt(dt_wide,
                id.vars       = "Sample",
                measure.vars  = ordre_final,
                variable.name = "Metric",
                value.name    = "Value")
dt_plot[, Metric := factor(Metric, levels = ordre_final)]

# ==========================================================================
# 4. Stats — moyenne % par métrique
# ==========================================================================
dt_stats <- dt_plot[dt_plot[Metric == "Brut", .(Sample, Val_Brut = Value)],
                    on = "Sample"]
dt_stats[, Pct := (Value / Val_Brut) * 100]
dt_stats <- dt_stats[, .(Mean_Pct = round(mean(Pct, na.rm = TRUE), 1)), by = Metric]

labels_dyn        <- setNames(
  paste0(dt_stats$Metric, " (", dt_stats$Mean_Pct, "%)"),
  as.character(dt_stats$Metric)
)
labels_dyn["Brut"] <- "Brut (100%)"

# ==========================================================================
# 5. Couleurs
# ==========================================================================
couleurs <- setNames(
  c("grey70", turbo(length(ordre_final) - 1L)),
  ordre_final
)

# Position numérique des échantillons pour les deux barres côte à côte
sample_levels <- unique(dt_plot$Sample)
dt_plot[, x_num := as.numeric(factor(Sample, levels = sample_levels))]

# ==========================================================================
# 6. Graphique -> PDF
# ==========================================================================
p <- ggplot(dt_plot) +
  geom_col(
    data  = dt_plot[Metric == "Brut"],
    aes(x = x_num - 0.2, y = Value, fill = Metric),
    width = 0.35
  ) +
  geom_col(
    data     = dt_plot[Metric != "Brut"],
    aes(x = x_num + 0.2, y = Value, fill = Metric),
    width    = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  scale_x_continuous(breaks = seq_along(sample_levels), labels = sample_levels) +
  scale_fill_manual(values = couleurs, breaks = ordre_final, labels = labels_dyn) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Sample", y = "Counts", fill = "Étapes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(PDF, width = 10, height = 7)
print(p)
dev.off()

message("✓ PDF     : ", PDF)
message("✓ Parquet : ", PARQUET)