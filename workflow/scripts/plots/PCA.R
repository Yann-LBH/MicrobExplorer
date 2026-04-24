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
# Configuration Snakemake
# ==========================================================================
try(snakemake, silent = TRUE)
if (!exists("snakemake")) {
  snakemake <- list(
    input   = list(
      parquet  = list.files("Statistics/5.Final_Intersection", "\\.parquet$", full.names = TRUE),
      metadata = "metadata/metadata.xlsx",
      physico  = "Physico-chimique/Thermophiles.xlsx"
    ),
    output  = list(
      pdf     = "Graphique/ACP/ACP_results.pdf",
      parquet = "Data/Parquet/pca_results.parquet",
      csv     = "Graphique/ACP/contributions.csv"
    ),
    params  = list(
      dim_x        = 1L,          # Dimension X (défaut 1)
      dim_y        = 2L,          # Dimension Y (défaut 2)
      physico_cols = c("pH", "AGV mg/l"),  # Colonnes physico à inclure
      n_top        = 10L          # Nombre de taxons top contributeurs
    )
  )
}

setDTthreads(0L)

dim_x        <- as.integer(snakemake$params$dim_x)
dim_y        <- as.integer(snakemake$params$dim_y)
physico_cols <- snakemake$params$physico_cols
n_top        <- as.integer(snakemake$params$n_top)
dim_names    <- paste0("Dim.", c(dim_x, dim_y))

dir.create(dirname(snakemake$output$pdf),     recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(snakemake$output$parquet), recursive = TRUE, showWarnings = FALSE)

# ==========================================================================
# 1. Chargement des données
# ==========================================================================

# Parquet : lecture multi-fichiers si nécessaire
files <- unlist(snakemake$input$parquet)
df_final <- rbindlist(lapply(files, read_parquet), use.names = TRUE, fill = TRUE)

# Metadata
meta <- as.data.table(read_excel(snakemake$input$metadata))
# colonnes attendues : sample_id, name, date

# Physico-chimique
df_physico <- as.data.table(read_excel(snakemake$input$physico))

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
                             .SDcols = physico_cols]

nb            <- estim_ncpPCA(df_physico_num, ncp.max = min(3L, ncol(df_physico_num) - 1L))
res_impute    <- imputePCA(df_physico_num, ncp = nb$ncp)
physico_complet <- as.data.frame(res_impute$completeObs)

# ==========================================================================
# 5. PCA
# ==========================================================================
final_tab_pca   <- cbind(donnees_clr, physico_complet, Digesteur = metadata$Digesteur)

n_taxons        <- ncol(donnees_clr)
n_physico       <- length(physico_cols)
idx_quanti_sup  <- (n_taxons + 1):(n_taxons + n_physico)
idx_quali_sup   <- n_taxons + n_physico + 1L

res_pca <- PCA(final_tab_pca,
               quanti.sup = idx_quanti_sup,
               quali.sup  = idx_quali_sup,
               graph      = FALSE)

# Top N taxons sur les dimensions choisies
contrib_sum    <- res_pca$var$contrib[, dim_x] + res_pca$var$contrib[, dim_y]
top_taxons     <- names(sort(contrib_sum, decreasing = TRUE)[seq_len(n_top)])
selection_finale <- c(top_taxons, physico_cols)

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
  axes        = c(dim_x, dim_y),       # dimensions configurables
  geom.ind    = "none",
  col.var     = "black",
  col.quanti.sup = "blue",
  select.var  = list(name = selection_finale),
  habillage   = idx_quali_sup,
  mean.point  = FALSE,
  repel       = TRUE,
  title       = sprintf("ACP CLR — Top %d taxons (Dim %d vs %d)", n_top, dim_x, dim_y)
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
pdf(snakemake$output$pdf, width = 12, height = 8)
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
coords_ind[, (physico_cols) := physico_complet]

write_parquet(coords_ind, snakemake$output$parquet)

# --- CSV contributions ---
fwrite(contrib_top, snakemake$output$csv, sep = ";")

message("✓ PDF       : ", snakemake$output$pdf)
message("✓ Parquet   : ", snakemake$output$parquet)
message("✓ Contrib   : ", snakemake$output$csv)