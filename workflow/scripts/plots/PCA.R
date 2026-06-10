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
library(FactoMineR)
library(factoextra)
library(compositions)
library(missMDA)
library(vegan)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Inputs
DATA        <- snakemake[["input"]][["data"]]
METADATA    <- snakemake[["input"]][["metadata"]]
PHYSICO     <- snakemake[["input"]][["physico"]]

# Outputs
PDF         <- snakemake[["output"]][["pdf"]]
PARQUET     <- snakemake[["output"]][["parquet"]]
CSV         <- snakemake[["output"]][["csv"]]

# Parameters
SHARED       <- snakemake[["params"]][["shared"]]
TOP_N        <- as.integer(snakemake[["params"]][["top_n"]])
POINT_SIZE   <- as.numeric(snakemake[["params"]][["point_size"]])
DIM_X        <- as.integer(snakemake[["params"]][["dim_x"]])
DIM_Y        <- as.integer(snakemake[["params"]][["dim_y"]])
PHYSICO_COLS <- snakemake[["params"]][["physico_cols"]]
#dim_names    <- paste0("Dim.", c(dim_x, dim_y))

# ==========================================================================
# 1. Data Import
# ==========================================================================

# Parquet input files
files <- unlist(DATA)
df_final <- rbindlist(lapply(files, function(f) read_parquet(f, col_select = c("Date", "Digesteur", "Contig_ID", "RPKM"))), use.names = TRUE, fill = TRUE)

# Metadata
meta <- as.data.table(read_excel(METADATA))
# colonnes attendues : sample_id, name, date

# Physico-chimique
df_physico <- as.data.table(read_excel(PHYSICO))

# ==========================================================================
# 2. Pivot wide + jointure metadata
# ==========================================================================
tableau_large <- dcast(
  df_final[, .(Date, Digesteur, Contig_ID, RPKM = as.numeric(RPKM))],
  Date + Digesteur ~ Contig_ID,
  value.var = "RPKM",
  fun.aggregate = sum,
  fill = 0
)

# Suppression colonne artefact si présente
tableau_large[, any_of("Other_NA") := NULL]

metadata <- tableau_large[, .(Date = as.factor(Date), Digesteur = as.factor(Digesteur))]
data_pca  <- tableau_large[, -(1:2)]

# ==========================================================================
# 3. Transformation CLR
# ==========================================================================
donnees_clr   <- as.data.frame(as.matrix(clr(data_pca + 1)))

# ==========================================================================
# 4. Physico-chimique : sélection dynamique + imputation
# ==========================================================================
# Colonnes configurables via params
df_physico_num <- df_physico[, lapply(.SD, function(x) as.numeric(gsub(",", ".", x))),
                             .SDcols = PHYSICO_COLS]

nb            <- estim_ncpPCA(df_physico_num, ncp.max = min(3L, ncol(df_physico_num) - 1L))
res_impute    <- imputePCA(df_physico_num, ncp = nb$ncp)
physico_complet <- as.data.frame(res_impute$completeObs)

# ==========================================================================
# 5. PCA
# ==========================================================================
final_tab_pca   <- cbind(donnees_clr, physico_complet, Digesteur = metadata$Digesteur)

n_taxons        <- ncol(donnees_clr)
n_physico       <- length(PHYSICO_COLS)
idx_quanti_sup  <- (n_taxons + 1):(n_taxons + n_physico)
idx_quali_sup   <- n_taxons + n_physico + 1L

res_pca <- PCA(final_tab_pca,
               quanti.sup = idx_quanti_sup,
               quali.sup  = idx_quali_sup,
               graph      = FALSE)

# Top N taxons sur les dimensions choisies
contrib_sum    <- res_pca$var$contrib[, DIM_X] + res_pca$var$contrib[, DIM_Y]
top_taxons     <- names(sort(contrib_sum, decreasing = TRUE)[seq_len(TOP_N)])
selection_finale <- c(top_taxons, PHYSICO_COLS)

# ==========================================================================
# 6. Graphiques
# ==========================================================================
n_dates     <- nlevels(metadata$Date)
n_digesteurs <- nlevels(metadata$Digesteur)

palette_digesteurs <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFC0CB")
formes_dates       <- c(16, 17, 15, 18, 8)

plot_variance <- fviz_eig(res_pca, addlabels = TRUE)

plot_acp <- fviz_pca_biplot(
  res_pca,
  axes        = c(DIM_X, DIM_Y),       # dimensions configurables
  geom.ind    = "none",
  col.var     = "black",
  col.quanti.sup = "blue",
  select.var  = list(name = selection_finale),
  habillage   = idx_quali_sup,
  mean.point  = FALSE,
  repel       = TRUE,
  title       = sprintf("ACP CLR — Top %d taxons (Dim %d vs %d)", TOP_N, DIM_X, DIM_Y)
) +
  geom_point(
    data = as.data.frame(res_pca$ind$coord),
    aes(x     = .data[[dim_names[1]]],
        y     = .data[[dim_names[2]]],
        shape = metadata$Date,
        color = metadata$Digesteur),
    size = 4, alpha = 0.8
  ) +
  scale_shape_manual(values = formes_dates[seq_len(n_dates)],
                     name   = "Dates",
                     labels = levels(metadata$Date)) +
  scale_color_manual(values = palette_digesteurs[seq_len(n_digesteurs)],
                     name   = "Digesteurs",
                     labels = levels(metadata$Digesteur)) +
  theme_minimal() +
  theme(legend.position = "right")

plot_acp$layers <- rev(plot_acp$layers)

# Contributions top N
contrib_top <- as.data.table(res_pca$var$contrib[top_taxons, ],
                             keep.rownames = "Taxon")

# ==========================================================================
# 7. Sorties
# ==========================================================================

# --- PDF compilé (ACP + Variance + Contributions) ---
pdf(PDF, width = 12, height = 8)
print(plot_variance)
print(plot_acp)
# Tableau contributions comme plot texte
grid::grid.newpage()
grid::grid.draw(gridExtra::tableGrob(contrib_top))
dev.off()

# --- Parquet pour Shiny ---
# Coordonnées individus + metadata + physico -> tout dans un seul fichier
coords_ind <- as.data.table(res_pca$ind$coord)
coords_ind[, Digesteur := metadata$Digesteur]
coords_ind[, Date      := metadata$Date]
coords_ind[, (PHYSICO_COLS) := physico_complet]

write_parquet(coords_ind, PARQUET)

# --- CSV contributions ---
fwrite(contrib_top, CSV, sep = ";")

message("✓ PDF       : ", PDF)
message("✓ Parquet   : ", PARQUET)
message("✓ Contrib   : ", CSV)