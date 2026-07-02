################################################################################
# Project : "MicrobExplorer"
# Script: "Plotting PCA Contigs with taxonomy"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

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

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Inputs
DATA <- as.character(snakemake@input[["data"]])
METADATA <- as.character(c(snakemake@input[["metadata"]])[1])
PHYSICO <- as.character(c(snakemake@input[["physico"]])[1])

# Outputs
PDF <- as.character(c(snakemake@output[["pdf"]])[1])
PARQUET <- as.character(c(snakemake@output[["parquet"]])[1])
XLSX <- as.character(c(snakemake@output[["xlsx"]])[1])

# Parameters
SHARED <- snakemake@params[["shared"]]
TOP_N <- as.integer(c(snakemake@params[["top_n"]])[1])
POINT_SIZE <- as.numeric(c(snakemake@params[["point_size"]])[1])
DIM_X <- as.integer(c(snakemake@params[["dim_x"]])[1])
DIM_Y <- as.integer(c(snakemake@params[["dim_y"]])[1])
PHYSICO_COLS <- as.character(snakemake@params[["physico_cols"]])

# Wildcards
SOURCE <- as.character(snakemake@wildcards[["source"]])

# --- Extraction et Conversion Type par Type ---

# 1. Les Textes (Character)
THEME   <- as.character(SHARED$theme %||% "theme_minimal")
PALETTE <- as.character(SHARED$palette %||% "Turbo")

# 2. Les Nombres Décimaux (Numeric / Double)
# Utile pour les dimensions des plots car on peut vouloir 8.5 ou 10.2
WIDTH  <- as.numeric(SHARED$width %||% 10)
HEIGHT <- as.numeric(SHARED$height %||% 8)

# 3. Les Tailles de Police / Nombres Entiers (Integer)
# Les tailles de texte sont généralement des entiers
TITLE_SIZE  <- as.integer(SHARED$title$size %||% 12)
LEGEND_SIZE <- as.integer(SHARED$legend$size %||% 10)
AXES_SIZE   <- as.integer(SHARED$axes$size %||% 10)

dim_names <- paste0("Dim.", c(DIM_X, DIM_Y))

# ==========================================================================
# 1. Data Import
# ==========================================================================

# Metadata loading
meta <- as.data.table(read_excel(METADATA))
setDT(meta)

meta[, sample_id := as.character(sample_id)]
meta[, name := as.character(name)]
meta[, date := as.character(date)]

all_sample_ids <- meta$sample_id

# Dynamically map column names based on source type
if (grepl("read", SOURCE, ignore.case = TRUE)) {
  ID_COL <- "read_id"
  ABUNDANCE_COL <- "cpm"
} else {
  ID_COL <- "contig_id"
  ABUNDANCE_COL <- "rpkm"
}

files <- unlist(DATA)

# Read all TSV files and map the sample_id using metadata matching
df_list <- lapply(files, function(f) {
  file_name_only <- basename(f)
  
  matched_id <- all_sample_ids[sapply(all_sample_ids, function(id) grepl(id, file_name_only, fixed = TRUE))]
  
  if (length(matched_id) == 0) {
    stop(paste("Could not map any sample_id from metadata to the file:", file_name_only))
  }
  
  dt <- fread(f, select = c(ID_COL, ABUNDANCE_COL))
  dt[, sample_id := matched_id[1]]
  dt[, (ABUNDANCE_COL) := as.numeric(get(ABUNDANCE_COL))]
  return(dt)
})

df_final <- rbindlist(df_list, use.names = TRUE, fill = TRUE)

# Merge date and name from metadata into df_final
df_final <- merge(df_final, meta[, .(sample_id, date, name)], by = "sample_id", all.x = TRUE)

# ==========================================================================
# 2. Pivot wide
# ==========================================================================
formula_dcast <- as.formula(paste("date + name ~", ID_COL))

tableau_large <- dcast(
  df_final,
  formula_dcast,
  value.var = ABUNDANCE_COL,
  fun.aggregate = sum,
  fill = 0
)

cols_to_remove <- intersect("Other_NA", names(tableau_large))
if (length(cols_to_remove) > 0) {
  tableau_large[, (cols_to_remove) := NULL]
}

# ==========================================================================
# 3. Safe Merge with Physico-chemical Data via date + name
# ==========================================================================

# Load physico Excel file
df_physico <- as.data.table(read_excel(PHYSICO))
setDT(df_physico)
df_physico[, date := as.character(date)]
df_physico[, name := as.character(name)]

# Clean and parse numeric variables from parameters sheet
df_physico_clean <- df_physico[, lapply(.SD, function(x) {
  clean_x <- gsub("[^0-9.,-]", "", as.character(x))
  clean_x <- gsub(",", ".", clean_x)
  return(as.numeric(clean_x))
}), .SDcols = PHYSICO_COLS]

# Add keys back to the cleaned physico dataset
df_physico_clean[, `:=`(date = df_physico$date, name = df_physico$name)]

# Align and merge physico data right into our wide data structure based on date and name
tableau_large <- merge(tableau_large, df_physico_clean, by = c("date", "name"), all.x = TRUE)

# Split metadata and numeric features
metadata <- tableau_large[, .(date = as.factor(date), name = as.factor(name))]

# Extract taxons data and physico data separately but aligned
taxon_cols <- setdiff(names(tableau_large), c("date", "name", PHYSICO_COLS))
data_pca <- tableau_large[, taxon_cols, with = FALSE]
physico_for_impute <- as.data.frame(tableau_large[, PHYSICO_COLS, with = FALSE])

# ==========================================================================
# 4. CLR Transformation & Imputation
# ==========================================================================
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
n_physico <- length(PHYSICO_COLS)
idx_quanti_sup <- (n_taxons + 1):(n_taxons + n_physico)
idx_quali_sup <- n_taxons + n_physico + 1L

res_pca <- PCA(final_tab_pca,
  quanti.sup = idx_quanti_sup,
  quali.sup  = idx_quali_sup,
  graph      = FALSE
)

contrib_sum <- res_pca$var$contrib[, DIM_X] + res_pca$var$contrib[, DIM_Y]
top_taxons <- names(sort(contrib_sum, decreasing = TRUE)[seq_len(TOP_N)])
selection_finale <- c(top_taxons, PHYSICO_COLS)

# ==========================================================================
# 6. Graphics Generation
# ==========================================================================
n_dates <- nlevels(metadata$date)
n_digesteurs <- nlevels(metadata$name)

palette_digesteurs <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFC0CB")
formes_dates <- c(16, 17, 15, 18, 8)

plot_variance <- fviz_eig(res_pca, addlabels = TRUE)

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
  title = sprintf("PCA CLR (%s) - Top %d features (Dim %d vs %d)", toupper(SOURCE), TOP_N, DIM_X, DIM_Y)
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
    values = formes_dates[seq_len(n_dates)],
    name = "Dates",
    labels = levels(metadata$date)
  ) +
  scale_color_manual(
    values = palette_digesteurs[seq_len(n_digesteurs)],
    name = "Digesteurs",
    labels = levels(metadata$name)
  ) +
  labs(x = paste("Dimension", DIM_X), y = paste("Dimension", DIM_Y)) +
  theme_minimal() +
  theme(legend.position = "right")

plot_acp$layers <- rev(plot_acp$layers)

contrib_top <- as.data.table(res_pca$var$contrib[top_taxons, ],
  keep.rownames = "Taxon"
)

# ==========================================================================
# 7. Outputs and File Generation
# ==========================================================================
pdf(PDF, width = 10, height = 8)
print(plot_variance)
print(plot_acp)

grid::grid.newpage()
grid::grid.draw(gridExtra::tableGrob(contrib_top))
dev.off()

coords_ind <- as.data.table(res_pca$ind$coord)
coords_ind[, name := metadata$name]
coords_ind[, date := metadata$date]
coords_ind[, (PHYSICO_COLS) := physico_complet]

write_parquet(coords_ind, PARQUET)

dir.create(dirname(XLSX), recursive = TRUE, showWarnings = FALSE)
write_xlsx(contrib_top, XLSX)

message("✓ PDF       : ", PDF)
message("✓ Parquet   : ", PARQUET)
message("✓ Excel     : ", XLSX)