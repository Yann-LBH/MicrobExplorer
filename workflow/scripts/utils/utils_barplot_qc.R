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
library(viridisLite)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Inputs
DATA <- snakemake[["input"]][["data"]]

# Outputs
PDF <- snakemake[["output"]][["pdf"]]
PARQUET <- snakemake[["output"]][["parquet"]]

# Parameters
ACTIVE_MODULES <- snakemake[["params"]][["active_modules"]]

# ==========================================================================
# 1. Loading and Dynamic Setup
# ==========================================================================
# English comment: Read pipeline summary dataset
dt_raw <- fread(DATA, sep = "|", strip.white = TRUE)

# ✅ FIX: Define 'ordre_final' dynamically based on active pipeline variables
# This handles the exact metrics order expected by the downsteam logic
ordre_final <- c("Brut", "P_Counted", "P_Trimmed", "P_RPKM", "P_RPKM_Filtered", "Intersection")

# ==========================================================================
# 2. In-place Loss Calculations
# ==========================================================================
dt_raw[, `:=`(
  P_Counted       = Brut - Counted,
  P_Trimmed       = Counted - Trimmed,
  P_RPKM          = Counted - RPKM,
  P_RPKM_Filtered = RPKM - RPKM_Filtered,
  Intersection    = RPKM_Filtered - Intersection
)]

# Intersect to make sure we only select columns that actually exist in the table
selected_cols <- c("Sample", intersect(ordre_final, names(dt_raw)))
dt_wide <- dt_raw[, .SD, .SDcols = selected_cols]

# ==========================================================================
# 3. Parquet Output (Wide format for Shiny app reference)
# ==========================================================================
write_parquet(dt_wide, PARQUET)

# Convert wide structure to long formatting for ggplot scaling
dt_plot <- melt(dt_wide,
  id.vars       = "Sample",
  measure.vars  = intersect(ordre_final, names(dt_wide)),
  variable.name = "Metric",
  value.name    = "Value"
)

# ✅ FIX: Filter out rows that do not contribute (Value == 0) to clean up the plot
dt_plot <- dt_plot[Value > 0]

# Enforce strict factor ordering based on our template vector
dt_plot[, Metric := factor(Metric, levels = ordre_final)]

# ==========================================================================
# 4. Statistics : Mean percentage per metric tracking
# ==========================================================================
# Extract baseline Brut value to normalize downstream percentages
dt_brut_baseline <- dt_plot[Metric == "Brut", .(Sample, Val_Brut = Value)]

dt_stats <- dt_plot[dt_brut_baseline, on = "Sample"]
dt_stats[, Pct := (Value / Val_Brut) * 100]
dt_stats <- dt_stats[, .(Mean_Pct = round(mean(Pct, na.rm = TRUE), 1)), by = Metric]

# Generate beautiful dynamic labels containing calculation metrics
labels_dyn <- setNames(
  paste0(dt_stats$Metric, " (", dt_stats$Mean_Pct, "%)"),
  as.character(dt_stats$Metric)
)
labels_dyn["Brut"] <- "Brut (100%)"

# ==========================================================================
# 5. Colors Mapping Configuration
# ==========================================================================
# ✅ FIX: turbo() is now safely executable because viridisLite is loaded
couleurs <- setNames(
  c("grey70", turbo(length(ordre_final) - 1L)),
  ordre_final
)

# Set dynamic X positions to split the raw vs losses bars side-by-side
sample_levels <- unique(dt_plot$Sample)
dt_plot[, x_num := as.numeric(factor(Sample, levels = sample_levels))]

# ==========================================================================
# 6. Graphics Generation -> PDF
# ==========================================================================
p <- ggplot(dt_plot) +
  geom_col(
    data  = dt_plot[Metric == "Brut"],
    aes(x = x_num - 0.2, y = Value, fill = Metric),
    width = 0.35
  ) +
  geom_col(
    data = dt_plot[Metric != "Brut"],
    aes(x = x_num + 0.2, y = Value, fill = Metric),
    width = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  scale_x_continuous(breaks = seq_along(sample_levels), labels = sample_levels) +
  scale_fill_manual(values = couleurs, breaks = ordre_final, labels = labels_dyn, na.value = "transparent") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Sample", y = "Counts", fill = "Steps") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Write output stream to PDF file target
pdf(PDF, width = 10, height = 7)
print(p)
dev.off()

message("✓ PDF     : ", PDF)
message("✓ Parquet : ", PARQUET)
