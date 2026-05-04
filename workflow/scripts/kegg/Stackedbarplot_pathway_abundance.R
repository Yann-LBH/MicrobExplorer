################################################################################
# Project : "MicrobExplorer"
# Script: "Plotting pathway abundance"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(data.table)
library(ggplot2)
library(arrow)
library(readxl) 

# --- Configuration ---
data_path      <- "Statistics/36_Kegg_merge_pathway_levels/"
metadata_path  <- "metadata/metadata.xlsx"
output_pdf     <- "Graphique/Barplot/StackedBarplots_by_Condition.pdf"
output_parquet <- "Data/Parquet/aggregated_data.parquet"

dir.create(dirname(output_pdf),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_parquet), recursive = TRUE, showWarnings = FALSE)

# Nombre de threads pour data.table (utilise tous les cœurs dispo)
setDTthreads(0L)

# --- 1. Chargement des données ---

meta <- as.data.table(read_excel(metadata_path))

# Lecture parallèle et agrégation immédiate avec rbindlist (pas de rbind successifs)
files <- list.files(data_path, pattern = "annotated_agreg_.*\\.tsv", full.names = TRUE)

all_data <- rbindlist(
  lapply(files, function(f) {
    dt <- fread(f,
                select      = c("level_3", "adj_standardization"),
                nThread     = getDTthreads(),  # parallélisme par fichier
                showProgress = FALSE)
    dt[, sample_id := gsub("annotated_agreg_|\\.tsv", "", basename(f), fixed = FALSE)]
    dt
  }),
  use.names = TRUE,
  fill      = TRUE      # robustesse si colonnes manquantes dans certains fichiers
)

# Jointure : meta à gauche pour garder uniquement les samples connus
setkey(all_data, sample_id)
setkey(meta,     sample_id)
all_data <- all_data[meta, nomatch = 0]   # inner join, all_data piloté par meta

write_parquet(all_data, output_parquet)

# --- 2. Agrégation (Top 10 + Others) ---
# Tout en data.table, zéro copie inutile

agg_data <- all_data[,
                     .(total_adj = sum(adj_standardization, na.rm = TRUE)),
                     by = .(name, date, level_3)
]

# Rank en place, sans colonne temporaire
agg_data[, category := {
  r <- frank(-total_adj, ties.method = "random")
  fifelse(r <= 10L, level_3, "Others")
}, by = .(name, date)]

# Supprime level_3 et total_adj intermédiaires dès que possible
final_plot_data <- agg_data[,
                            .(sum_adj = sum(total_adj)),
                            by = .(name, date, category)
]
rm(agg_data, all_data)   # libère la RAM immédiatement

final_plot_data[,
                percentage := (sum_adj / sum(sum_adj)) * 100,
                by = .(name, date)
]

# --- 3. Facteurs & couleurs ---
top_pathways <- sort(unique(final_plot_data$category[final_plot_data$category != "Others"]))
levels_order <- c(top_pathways, "Others")
final_plot_data[, category := factor(category, levels = levels_order)]

# Palette turbo sans dépendance à scales::
# grDevices::hcl.colors est base R, mais turbo vient de viridisLite (dépendance de ggplot2)
# -> viridisLite est donc déjà chargé via ggplot2, pas de lib supplémentaire
my_colors        <- viridisLite::viridis(length(top_pathways), option = "turbo")
my_colors        <- c(my_colors, "#000000")
names(my_colors) <- levels_order

# --- 4. Génération PDF ---
# Pré-split par condition : évite de filtrer le DT à chaque itération
split_data <- split(final_plot_data, by = "name", keep.by = TRUE)

pdf(output_pdf, width = 12, height = 8)

lapply(split_data, function(df_cond) {
  p <- ggplot(df_cond, aes(x = as.factor(date), y = percentage, fill = category)) +
    geom_bar(stat = "identity",
             position  = position_stack(reverse = TRUE),
             color     = "white",
             linewidth = 0.05) +
    scale_fill_manual(values = my_colors) +
    labs(
      title    = paste("Condition:", df_cond$name[1L]),
      subtitle = "Top 10 abundance pathways",
      x        = "Date",
      y        = "Relative Abundance (%)",
      fill     = "Pathways"
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      legend.text  = element_text(size = 7),
      panel.grid.major.x = element_blank()
    )
  print(p)
})

dev.off()
message("Task complete. PDF: ", output_pdf)