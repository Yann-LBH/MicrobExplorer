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
DATA <- as.character(snakemake@input[["data"]])

# Outputs
PDF <- as.character(snakemake@output[["pdf"]])[1]
PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Parameters
STEPS_CONFIG <- snakemake@params[["steps_config"]]
ACTIVE_MODULES <- as.character(snakemake@params$active_modules)

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]
# ==========================================================================
# 1. Loading and Dynamic Setup
# ==========================================================================
# Read pipeline summary dataset from Parquet binary format
dt_raw <- as.data.table(read_parquet(DATA))

base_col <- if ("extracted" %in% names(dt_raw)) "extracted" else "brut"

if ("Sample" %in% names(dt_raw)) setnames(dt_raw, "Sample", "sample")

chain <- names(STEPS_CONFIG)

if (is.null(chain) || length(chain) < 2) {
  stop(sprintf("❌ Error: steps_config invalide ou vide pour source [%s].", SOURCE))
}
# ==========================================================================
# 2. In-place Loss Calculations
# ==========================================================================
ordre_final <- character(0)

for (i in seq_len(length(chain) - 1)) {
  prev_step <- chain[i]
  curr_step <- chain[i + 1]
  loss_col  <- paste0("p_", curr_step)

  if (prev_step %in% names(dt_raw) && curr_step %in% names(dt_raw)) {
    dt_raw[, (loss_col) := get(prev_step) - get(curr_step)]
    ordre_final <- c(ordre_final, curr_step, loss_col)
  } else {
    missing <- setdiff(c(prev_step, curr_step), names(dt_raw))
    warning(sprintf("⚠️ [%s] Étape ignorée (%s → %s) : colonne(s) manquante(s) : %s",
                     SOURCE, prev_step, curr_step, paste(missing, collapse = ", ")))
  }
}

ordre_final <- c(base_col, ordre_final)
# ==========================================================================
# 3. Parquet Output & Hard-Secured Reshaping
# ==========================================================================
available_metrics <- intersect(ordre_final, names(dt_raw))
if (length(available_metrics) == 0 && "brut" %in% names(dt_raw)) {
  available_metrics <- "brut"
}

selected_cols <- c("sample", available_metrics)
dt_wide <- dt_raw[, .SD, .SDcols = selected_cols]

# Save the wide structure to output target
write_parquet(dt_wide, PARQUET)

# Standard base R dataframe casting to ensure stable pivoting environment
df_standard <- as.data.frame(dt_wide)
df_long <- reshape(
  df_standard, 
  varying       = available_metrics, 
  v.names       = "value",
  timevar       = "variable",
  times         = available_metrics,
  direction     = "long",
  new.row.names = NULL
)

dt_plot <- as.data.table(df_long)
if ("id" %in% names(dt_plot)) dt_plot[, id := NULL]

setnames(dt_plot, "sample", "Sample", skip_absent = TRUE)
dt_plot[, variable := as.character(variable)]
dt_plot[, value    := as.numeric(value)]

# ==========================================================================
# 4. Debugging Condition for Non-Impact Steps
# ==========================================================================
if (nrow(dt_plot) > 0) {
  metric_impacts <- dt_plot[, .(Total_Loss = sum(value, na.rm = TRUE)), by = variable]
  # Protect the dynamic baseline instead of hardcoded "brut"
  zero_impact_steps <- metric_impacts[Total_Loss == 0 & variable != base_col, variable]
  dt_plot <- dt_plot[!variable %in% zero_impact_steps]
  
  # Crucial fix: Ensure base_col is NEVER filtered out even if value is 0 or NA
  dt_plot <- dt_plot[value > 0 | variable == base_col | is.na(value)]
}

# Clean up secondary baseline data if we are tracking an upgraded unit scale (e.g. genes)
if (base_col == "extracted") {
  dt_plot <- dt_plot[variable != "brut"]
}

active_order <- intersect(ordre_final, unique(dt_plot$variable))
if (length(active_order) > 0) {
  dt_plot[, variable := factor(variable, levels = active_order)]
}

# ==========================================================================
# 5. Statistics & Labels Tuning
# ==========================================================================
# Dynamic baseline matching using our global parameter
dt_brut_baseline <- dt_plot[variable == base_col, .(Sample, Val_Brut = value)]

# Initialize labels_dyn with default names for ALL active metrics
labels_dyn <- setNames(as.character(active_order), active_order)

if (nrow(dt_plot[variable != base_col]) > 0 && nrow(dt_brut_baseline) > 0) {
  dt_stats <- dt_plot[dt_brut_baseline, on = "Sample"]
  dt_stats[, Pct := ifelse(Val_Brut > 0, (value / Val_Brut) * 100, 0)]
  dt_stats <- dt_stats[, .(Mean_Pct = round(mean(Pct, na.rm = TRUE), 1)), by = variable]
  
  # Update dynamic labels only for metrics that have calculated stats
  for (i in seq_len(nrow(dt_stats))) {
    m <- as.character(dt_stats$variable[i])
    labels_dyn[m] <- paste0(m, " (", dt_stats$Mean_Pct[i], "%)")
  }
}

# Enforce explicit label for the baseline
labels_dyn[base_col] <- paste0(base_col, " (100%)")
labels_dyn <- labels_dyn[active_order]

# ==========================================================================
# 6. Colors Mapping Configuration
# ==========================================================================
# Safety check to avoid palette generation errors if there are too few variables
num_colors <- max(1L, length(active_order) - 1L)
if (length(active_order) > 0) {
  couleurs <- setNames(
    c("grey70", turbo(num_colors))[1:length(active_order)],
    active_order
  )
} else {
  couleurs <- character(0)
}

sample_levels <- unique(dt_plot$Sample)
dt_plot[, x_num := as.numeric(factor(Sample, levels = sample_levels))]

# ==========================================================================
# 7. Graphics Generation -> PDF
# ==========================================================================
p <- ggplot(dt_plot) +
  geom_col(
    data  = dt_plot[variable == base_col],
    aes(x = x_num - 0.2, y = value, fill = variable),
    width = 0.35
  ) +
  geom_col(
    data     = dt_plot[variable != base_col],
    aes(x = x_num + 0.2, y = value, fill = variable),
    width    = 0.35,
    position = position_stack(reverse = TRUE)
  ) +
  scale_x_continuous(breaks = seq_along(sample_levels), labels = sample_levels) +
  scale_fill_manual(values = couleurs, breaks = active_order, labels = labels_dyn, na.value = "transparent") +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Sample", y = "Counts", fill = "Steps") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(PDF, width = 10, height = 7)
print(p)
dev.off()

message("✓ PDF     : ", PDF)
message("✓ Parquet : ", PARQUET)