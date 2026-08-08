################################################################################
# Project : "MicrobExplorer"
# Script: "Plotting PCA Contigs with taxonomy"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

suppressPackageStartupMessages({
# Libraries CRAN
  library(data.table)
  library(ggplot2)
  library(arrow)
  library(readxl)
  library(writexl)
  library(rlang)
  library(FactoMineR)
  library(factoextra)
  library(compositions)
  library(missMDA)
  library(vegan)
})

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================
source("workflow/scripts/utils/utils_title_resolver.R")
source("workflow/scripts/utils/utils_io.R")

# Inputs
DATA <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])[1]
PHYSICO <- as.character(snakemake@input[["physico"]])[1]
TAXONOMY <- as.character(snakemake@input[["taxonomy"]])[1]

# Outputs
PDF <- as.character(snakemake@output[["pdf"]])[1]
PARQUET <- as.character(snakemake@output[["parquet"]])[1]
XLSX <- as.character(snakemake@output[["xlsx"]])[1]

# Shared plots features
SHARED      <- snakemake@params[["shared"]]
THEME       <- as.character(SHARED$theme) %||% "theme_minimal"
PDF_SIZE    <- as.numeric(SHARED$pdf_size) %||% c(10, 8)
TITLE_SIZE  <- as.integer(SHARED$title_size) %||% 12
SUBTITLE_SIZE <- as.integer(SHARED$subtitle_size) %||% 10
LEGEND_SIZE <- as.integer(SHARED$legend_size) %||% 10
AXES_SIZE   <- as.integer(SHARED$axes_size) %||% 10
RANK        <- tolower(as.character(snakemake@params[["rank"]])[1])

# Parameters
TOP_N             <- as.integer(snakemake@params[["top_n"]])[1]
POINT_SIZE        <- as.numeric(snakemake@params[["point_size"]])[1]
DIM_X             <- as.integer(snakemake@params[["dim_x"]])[1]
DIM_Y             <- as.integer(snakemake@params[["dim_y"]])[1]
PHYSICO_COL      <- as.character(snakemake@params[["physico_col"]])
STAND_COL         <- as.character(snakemake@params[["stand_col"]])[1]

dim_names <- paste0("Dim.", c(DIM_X, DIM_Y))

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

RESOLVED_TITLE    <- resolve_text("TITLE_PCA", source = toupper(SOURCE), top_n = TOP_N)
RESOLVED_SUBTITLE <- resolve_text("SUBTITLE_PCA", dim_x = DIM_X, dim_y = DIM_Y, rank = RANK)

# ==========================================================================
# 1. Data Import & Redirection automatique KEGG
# ==========================================================================
if (is.na(RANK) || RANK == "" || RANK == "null") {
  stop("❌ Erreur critique : Le paramètre 'rank' est obligatoire dans la configuration Snakemake et ne peut pas être vide.")
}

# Automatic Detection and Redirection for KEGG
is_kegg <- grepl("kegg", SOURCE, ignore.case = TRUE)

if (is_kegg) {
  target_id <- "kegg_id"
} else {
  if (grepl("reads", SOURCE, ignore.case = TRUE)) {
    target_id <- "read_id" 
  } else if (grepl("contigs", SOURCE, ignore.case = TRUE)) {
    target_id <- "contig_id"
  } else {
    stop(sprintf("❌ Error: Unknown SOURCE value [%s]. Expected 'kegg', 'reads', or 'contigs'.", SOURCE))
  }
}

meta_dt <- load_metadata(METADATA)

all_data <- load_tsv_dir_dynamic(
  paths       = DATA,
  meta_dt     = meta_dt,
  select_cols = unique(c(target_id, "sample_id", RANK, STAND_COL))
)

# ==========================================================================
# TAXONOMY loading and merging
# ==========================================================================

dt_taxo <- fread(TAXONOMY)

# 3. Conditional merge with deduplication using data.table syntax
if (SOURCE == "reads") {
  dt_taxo_clean <- unique(dt_taxo, by = "tax_id")
  dt_merged     <- merge(all_data, dt_taxo_clean, by.x = "read_id", by.y = "tax_id", all.x = TRUE)
  tax_ranks     <- c("tax_id", "scientific_name", "domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")

} else if (SOURCE == "contigs") {
  dt_taxo_clean <- unique(dt_taxo, by = "contig_id")
  dt_merged     <- merge(all_data, dt_taxo_clean, by.x = "contig_id", by.y = "contig_id", all.x = TRUE)
  tax_ranks     <- c("domain", "phylum", "class", "order", "family", "genus", "species")

} else { # KEGG or other sources
  dt_taxo_clean <- unique(dt_taxo, by = "ko")
  dt_merged     <- merge(all_data, dt_taxo_clean, by.x = "kegg_id", by.y = "ko", all.x = TRUE)
  tax_ranks     <- c("ec_number", "level_1", "level_2", "level_3", "gene_description")
}

# 4. Clean up unassigned taxonomy in-place (replaces mutate + across)
tax_ranks_present <- intersect(tax_ranks, names(dt_merged))
for (col in tax_ranks_present) {
  set(dt_merged, i = which(is.na(dt_merged[[col]])), j = col, value = "Unclassified")
}

# ==========================================================================
# Label Processing and Wide Rotation
# ==========================================================================
# Format IDs as character and clean missing values
dt_merged[, (target_id) := as.character(get(target_id))]
dt_merged <- dt_merged[!is.na(get(target_id)) & get(target_id) != "" & get(target_id) != "NA"]

if (is_kegg) {
  # Cleaning the description column
  dt_merged[is.na(get(RANK)) | get(RANK) == "" | get(RANK) == "Unassigned", (RANK) := "Unknown Function"]
  # Kegg Label Format : "K00163 | Pyruvate dehydrogenase"
  feature_col <- "combined_label"
  dt_merged[, (feature_col) := paste(get(target_id), get(RANK), sep = " | ")]
} else {
  dt_merged[is.na(get(RANK)) | get(RANK) == "", (RANK) := get(target_id)]
  feature_col <- RANK
}

# Merge date and name from metadata into df_final
dt_merged <- merge(dt_merged, meta_dt[, .(sample_id, date, name)], by = "sample_id", all.x = TRUE)

formula_dcast <- as.formula(paste("date + name ~", feature_col))

tableau_large <- dcast(
  dt_merged,
  formula_dcast,
  value.var = STAND_COL,
  fun.aggregate = function(x) sum(x, na.rm = TRUE),
  fill = 0
)

# SECURITY
n_na <- sum(sapply(tableau_large, function(col) sum(is.na(col))))
if (n_na > 0) message(sprintf("Replacing %d remaining NA(s) with 0 after dcast.", n_na))
for (col in names(tableau_large)) {
  set(tableau_large, which(is.na(tableau_large[[col]])), col, 0)
}

# ==========================================================================
# 3. Merger of Physical Chemistry & NA Allocation
# ==========================================================================

# Load physico Excel file
dt_physico <- as.data.table(read_excel(PHYSICO))
dt_physico[, date := as.character(date)]
dt_physico[, name := as.character(name)]

# Clean and parse numeric variables from parameters sheet
dt_physico_clean <- dt_physico[, lapply(.SD, function(x) {
  clean_x <- gsub("[^0-9.,-]", "", as.character(x))
  clean_x <- gsub(",", ".", clean_x)
  result <- as.numeric(clean_x)
  # Security if number format is "1,235.5"
  n_new_na <- sum(is.na(result) & !is.na(x) & trimws(as.character(x)) != "")
  if (n_new_na > 0) {
    warning(sprintf("%d value(s) could not be parsed as numeric and became NA.", n_new_na))
  }
  result
}), .SDcols = PHYSICO_COL]

# Add keys back to the cleaned physico dataset
dt_physico_clean[, `:=`(date = dt_physico$date, name = dt_physico$name)]

# Align and merge physico data right into our wide data structure based on date and name
tableau_large <- merge(tableau_large, dt_physico_clean, by = c("date", "name"), all.x = TRUE)

# Split metadata and numeric features
metadata <- tableau_large[, .(date = as.factor(date), name = as.factor(name))]

# Extract taxons data and physico data separately but aligned
taxon_cols <- setdiff(names(tableau_large), c("date", "name", PHYSICO_COL))
data_pca <- tableau_large[, taxon_cols, with = FALSE]
physico_for_impute <- as.data.frame(tableau_large[, PHYSICO_COL, with = FALSE])

# ==========================================================================
# 4. CLR Transformation & Imputation
# ==========================================================================
# Security replace NA by 0
data_pca[is.na(data_pca)] <- 0
donnees_clr <- as.data.frame(as.matrix(clr(data_pca + 1)))

# Dynamic calculation for safe missMDA imputation bounds
max_possible_ncp <- min(nrow(physico_for_impute) - 2L, ncol(physico_for_impute) - 1L)

if (max_possible_ncp < 1L) {
  chosen_ncp <- 0L
} else {
  ncp_bound <- min(2L, max_possible_ncp)
  nb <- estim_ncpPCA(physico_for_impute, ncp.max = ncp_bound)
  chosen_ncp <- nb$ncp
}

res_impute <- imputePCA(physico_for_impute, ncp = chosen_ncp)
physico_complet <- as.data.frame(res_impute$completeObs)

# ==========================================================================
# 5. PCA Execution
# ==========================================================================
final_tab_pca <- cbind(donnees_clr, physico_complet, name = metadata$name)

n_taxons <- ncol(donnees_clr)
n_physico <- length(PHYSICO_COL)
idx_quanti_sup <- (n_taxons + 1):(n_taxons + n_physico)
idx_quali_sup <- n_taxons + n_physico + 1L

res_pca <- PCA(final_tab_pca,
  scale.unit = FALSE,
  ncp = max(DIM_X, DIM_Y, 5), # Security for dimensions
  quanti.sup = idx_quanti_sup,
  quali.sup  = idx_quali_sup,
  graph      = FALSE
)

contrib_sum <- res_pca$var$contrib[, DIM_X] + res_pca$var$contrib[, DIM_Y]
actual_top_n <- min(TOP_N, length(contrib_sum))
top_taxons <- names(sort(contrib_sum, decreasing = TRUE)[seq_len(actual_top_n)])
selection_finale <- c(top_taxons, PHYSICO_COL)

# ==========================================================================
# 6. Graphics Generation & Annotation
# ==========================================================================
n_dates <- nlevels(metadata$date)
n_samples <- nlevels(metadata$name)

color_blind_friendly_base <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#F0E442", # Yellow
  "#0072B2", # Blue
  "#D55E00", # Vermilion
  "#CC79A7", # Reddish Purple
  "#000000"  # Black
)

if (n_samples <= length(color_blind_friendly_base)) {
  palette_samples <- color_blind_friendly_base[seq_len(n_samples)]
} else {
  palette_samples <- grDevices::hcl.colors(n_samples, palette = "Qualitative")
}

base_shapes <- c(16, 17, 15, 18, 8, 1, 2, 0, 5, 6, 9, 10, 12, 13, 14)
if (n_dates <= length(base_shapes)) {
  shapes_dates <- base_shapes[seq_len(n_dates)]
} else {
  shapes_dates <- rep(base_shapes, length.out = n_dates)
}

plot_variance <- fviz_eig(res_pca, addlabels = TRUE) + 
  get(THEME)() +
  theme(plot.title = element_text(size = TITLE_SIZE))

# Defining dimensions for ggplot
dim_names <- paste0("Dim.", c(DIM_X, DIM_Y))

plot_acp <- fviz_pca_biplot(
  res_pca,
  axes = c(DIM_X, DIM_Y),
  geom.ind = "none",
  col.var = "black",
  col.quanti.sup = "blue",
  select.var = list(name = selection_finale),
  habillage = idx_quali_sup,
  mean.point = FALSE,
  repel = TRUE,
  title = RESOLVED_TITLE,
  subtitle = RESOLVED_SUBTITLE
) +
  geom_point(
    data = as.data.frame(res_pca$ind$coord),
    aes(
      x = .data[[dim_names[1]]],
      y = .data[[dim_names[2]]],
      shape = metadata$date,
      color = metadata$name
    ),
    size = POINT_SIZE, alpha = 0.8
  ) +
  scale_shape_manual(
    values = shapes_dates[seq_len(n_dates)],
    name = "Dates",
    labels = levels(metadata$date)
  ) +
  scale_color_manual(
    values = palette_samples[seq_len(n_samples)],
    name = "Samples",
    labels = levels(metadata$name)
  ) +
  labs(x = paste("Dimension", DIM_X), y = paste("Dimension", DIM_Y)) +
  get(THEME)() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = TITLE_SIZE, face = "bold"),
    plot.subtitle = element_text(size = SUBTITLE_SIZE),
    legend.text = element_text(size = LEGEND_SIZE),
  )

plot_acp$layers <- rev(plot_acp$layers)

contrib_top <- as.data.table(res_pca$var$contrib[top_taxons, ],
  keep.rownames = "Taxon"
)

# ==========================================================================
# 7. Outputs and File Generation
# ==========================================================================
pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])
on.exit(if (names(dev.cur()) != "null device") dev.off())
print(plot_variance)
print(plot_acp)

grid::grid.newpage()
grid::grid.draw(gridExtra::tableGrob(contrib_top))

coords_ind <- as.data.table(res_pca$ind$coord)
coords_ind[, name := metadata$name]
coords_ind[, date := metadata$date]
coords_ind[, (PHYSICO_COL) := physico_complet]

write_parquet(coords_ind, PARQUET)

write_xlsx(contrib_top, XLSX)

message("✓ PDF       : ", PDF)
message("✓ Parquet   : ", PARQUET)
message("✓ Excel     : ", XLSX)