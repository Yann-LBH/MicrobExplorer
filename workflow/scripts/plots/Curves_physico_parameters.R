# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Curves_phisyco_parameters.R
# PURPOSE : Processing and Multi-page PDF Plotting of Physicochemical Data
# AUTHOR  : Yann Le Bihan
# DATE    : 2026-09-03
# LINK    : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. METADATA & LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  # Libraries CRAN
  library(readxl)
  library(data.table)
  library(purrr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(arrow)
})

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# Helper: Clean and Convert Physicochemical Columns
cleaning_dataframe <- function(data, physicochemical_parameter) {
  df_clean <- data %>%
    mutate(across(
      all_of(physicochemical_parameter),
      ~ as.numeric(gsub(",", ".", as.character(.x)))
    )) %>%
    mutate(date = as.Date(date))
  return(df_clean)
}

# Helper: Generate Individual ggplot Curves with Ribbon
plot_generation_loop <- function(data, var_y, var_moy, var_sd, label_y,
                                 y_limites = NULL, ligne_seuil = NULL) {
  # Dataset for ribbon (variability)
  data_ribbon <- data %>%
    select(date, all_of(c(var_moy, var_sd))) %>%
    distinct() %>%
    arrange(date)

  remplissage_stats <- c("Ecart-Type" = "#8491B4")
  guide_style       <- list(title.position = "top", title.hjust = 0.5)

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
      colour   = do.call(guide_legend, c(list(order = 1), guide_style)),
      fill     = do.call(guide_legend, c(list(override.aes = list(alpha = 0.3), order = 2), guide_style)),
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

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_PHYSICO <- as.character(snakemake@input[["physico_params"]])[1]

# Outputs
OUT_PDF     <- as.character(snakemake@output[["pdf"]])[1]
OUT_PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Layout parameters
plots_per_page <- 6
nb_cols        <- 2
nb_rows        <- 3

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
if (is.null(IN_PHYSICO) || !file.exists(IN_PHYSICO)) {
  stop(sprintf("❌ Critical Error: Physicochemical input file '%s' does not exist.", IN_PHYSICO))
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading raw physicochemical Excel data...")
data_raw <- read_excel(IN_PHYSICO)

fixed_cols <- c("date", "name", "condition", "commentaire")
physicochemical_parameter <- setdiff(names(data_raw), fixed_cols)

if (length(physicochemical_parameter) == 0) {
  stop("❌ Critical Error: No physicochemical parameter columns found in input file.")
}

message(sprintf("✓ File loaded successfully: %d rows, %d parameters identified.", 
                nrow(data_raw), length(physicochemical_parameter)))

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PROCESSING FUNCTIONS
# ------------------------------------------------------------------------------
message("INFO: Cleaning data and calculating summary statistics...")
df_clean <- cleaning_dataframe(data_raw, physicochemical_parameter)

# Calcul des moyennes et écarts-types par date
df_stats <- df_clean %>%
  group_by(date) %>%
  summarise(across(all_of(physicochemical_parameter),
    list(
      Moyenne   = ~ mean(., na.rm = TRUE),
      EcartType = ~ sd(., na.rm = TRUE)
    ),
    .names = "{col}_{fn}"
  ), .groups = "drop")

df_final <- df_clean %>% left_join(df_stats, by = "date")

# ------------------------------------------------------------------------------
# 6. EXECUTION CORE & GRAPHICS GENERATION
# ------------------------------------------------------------------------------
message("INFO: Generating multipage PDF report...")
all_conditions <- na.omit(unique(df_final$condition))

pdf(PDF, width = 11, height = 8.5)
on.exit(if (names(dev.cur()) != "null device") dev.off())

walk(all_conditions, function(cond) {
  data_subset <- df_final %>% filter(condition == cond)

  # Generate all plots for this condition
  plot_list <- map(physicochemical_parameter, function(col) {
    if (sum(!is.na(data_subset[[col]])) > 0) {
      plot_generation_loop(
        data    = data_subset,
        var_y   = col,
        var_moy = paste0(col, "_Moyenne"),
        var_sd  = paste0(col, "_EcartType"),
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
          title    = paste("Analyse Condition :", cond),
          subtitle = paste("Page", i, "sur", length(pages), "|", length(plot_list), "paramètres"),
          theme    = theme(plot.title = element_text(size = 18, face = "bold"))
        ) +
        plot_layout(guides = "collect") &
        theme(
          legend.position = "bottom",
          legend.text     = element_text(size = 10),
          axis.title      = element_text(size = 10, face = "bold"),
          axis.text       = element_text(size = 9)
        )

      print(combined_plot)
    })
  }
})

# Fallback if no plot has been generated
if (length(combined_plot) == 0) {
  plot.new()
  msg <- "No plots has been drawn"
  text(0.5, 0.5, msg, cex = 1.1)
}

if (names(dev.cur()) != "null device") dev.off()

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (length(combined_plot) > 0) {
  arrow::write_parquet(df_final, OUT_PARQUET)
  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: Empty output: No plots has been drawn")
}
