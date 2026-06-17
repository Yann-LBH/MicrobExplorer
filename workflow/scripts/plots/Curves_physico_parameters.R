################################################################################
# Project : "MicrobExplorer"
# Script: "Plotting physico-chemical parameters"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
library(readxl)
library(data.table)
library(purrr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(arrow)

# ==========================================================================
# Configuration Snakemake
# ==========================================================================

# Inputs
# MODIFICATION : Sécurisation absolue du type pour read_excel
PHYSICO <- as.vector(as.character(snakemake@input$physico_params))

# Outputs
PDF <- as.character(snakemake@output$pdf)
PARQUET <- as.character(snakemake@output$parquet)

# ==========================================================================
# Functions
# ==========================================================================

cleaning_dataframe <- function(data, physicochemical_parameter) {
  df_clean <- data %>%
    mutate(across(
      all_of(physicochemical_parameter),
      ~ as.numeric(gsub(",", ".", as.character(.x)))
    )) %>%
    mutate(date = as.Date(date))
  return(df_clean)
}

plot_generation_loop <- function(data, var_y, var_moy, var_sd, label_y,
                                 y_limites = NULL, ligne_seuil = NULL) {
  # Dataset for ribbon (variability)
  data_ribbon <- data %>%
    select(date, all_of(c(var_moy, var_sd))) %>%
    distinct() %>%
    arrange(date)

  remplissage_stats <- c("Ecart-Type" = "#8491B4")
  guide_style <- list(title.position = "top", title.hjust = 0.5)

  p <- ggplot(data, aes(x = date, y = .data[[var_y]], group = name, color = name)) +
    # Ribbon
    geom_ribbon(
      data = data_ribbon, aes(
        x = date, ymin = .data[[var_moy]] - .data[[var_sd]],
        ymax = .data[[var_moy]] + .data[[var_sd]],
        fill = "Ecart-Type", group = 1
      ), inherit.aes = FALSE,
      alpha = 0.3, colour = NA
    ) +
    # Lines
    geom_line(linewidth = 1.2) +
    geom_line(aes(y = .data[[var_moy]], linetype = "Moyenne"), color = "Black", linewidth = 1.7) +
    # Scales
    scale_colour_viridis_d(name = "Digesteurs", option = "Turbo") +
    scale_fill_manual(name = "Variabilité", values = remplissage_stats) +
    scale_linetype_manual(name = "Tendance", values = c("Moyenne" = "solid")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    # Legend style
    guides(
      colour = do.call(guide_legend, c(list(order = 1), guide_style)),
      fill = do.call(guide_legend, c(list(override.aes = list(alpha = 0.3), order = 2), guide_style)),
      linetype = do.call(guide_legend, c(list(order = 3), guide_style))
    ) +
    labs(x = "Date", y = label_y) +
    theme_light()

  # Threshold line
  if (!is.null(ligne_seuil)) {
    p <- p + geom_hline(yintercept = ligne_seuil, color = "red", linewidth = 1, linetype = "dashed") +
      annotate("text",
        x = min(data$date), y = ligne_seuil,
        label = "Seuil", color = "red", vjust = -0.5, size = 3.5
      )
  }

  if (!is.null(y_limites)) {
    p <- p + coord_cartesian(ylim = y_limites)
  }

  return(p)
}

# ==========================================================================
# Data Processing & Parquet Export
# ==========================================================================

# 1. Load data
data_raw <- read_excel(PHYSICO)
fixed_cols <- c("date", "name", "condition", "commentaire")
physicochemical_parameter <- setdiff(names(data_raw), fixed_cols)

# 2. Clean and Compute Stats
df_clean <- cleaning_dataframe(data_raw, physicochemical_parameter)

df_stats <- df_clean %>%
  group_by(date) %>%
  summarise(across(all_of(physicochemical_parameter),
    list(
      Moyenne = ~ mean(., na.rm = TRUE),
      EcartType = ~ sd(., na.rm = TRUE)
    ),
    .names = "{col}_{fn}"
  ), .groups = "drop")

df_final <- df_clean %>% left_join(df_stats, by = "date")

# 3. Export to Parquet
write_parquet(df_final, PARQUET)
message("\n✓ PARQUET : Curves_physico_parameters Successfully generated.")

# ==========================================================================
# Optimized Plotting Loop
# ==========================================================================

# Layout settings
plots_per_page <- 6
nb_cols <- 2
nb_rows <- 3

all_conditions <- na.omit(unique(df_final$condition))

# MODIFICATION CRITIQUE : On ouvre le fichier PDF unique ICI avant de commencer à boucler
pdf(PDF, width = 11, height = 8.5)

walk(all_conditions, function(cond) {
  data_subset <- df_final %>% filter(condition == cond)

  # Generate all plots for this condition
  plot_list <- map(physicochemical_parameter, function(col) {
    if (sum(!is.na(data_subset[[col]])) > 0) {
      plot_generation_loop(
        data = data_subset,
        var_y = col,
        var_moy = paste0(col, "_Moyenne"),
        var_sd = paste0(col, "_EcartType"),
        label_y = col
      )
    } else {
      NULL
    }
  })

  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) > 0) {
    # Split into pages
    pages <- split(plot_list, ceiling(seq_along(plot_list) / plots_per_page))

    walk2(pages, seq_along(pages), function(page_plots, i) {
      combined_plot <- wrap_plots(page_plots, ncol = nb_cols, nrow = nb_rows) +
        plot_annotation(
          title = paste("Analyse Condition :", cond),
          subtitle = paste("Page", i, "sur", length(pages), "|", length(plot_list), "paramètres"),
          theme = theme(plot.title = element_text(size = 18, face = "bold"))
        ) +
        plot_layout(guides = "collect") &
        theme(
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 9)
        )

      print(combined_plot)
    })
  }
})
dev.off()

message("\n✓ PDF : Curves_physico_parameters Successfully generated.", PDF)
