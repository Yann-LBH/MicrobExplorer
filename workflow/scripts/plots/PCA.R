# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : PCA.R
# PURPOSE : CLR Transformation, Imputation & PCA Biplot integrating Phyloseq 
#           & Physico-Chemistry
# AUTHOR  : Yann Le Bihan
# DATE    : 2026-09-03
# LINK    : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. METADATA & LIBRARIES
# ------------------------------------------------------------------------------
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

source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_PHYLOSEQ <- as.character(snakemake@input[["phyloseq_obj"]])[1]
IN_PHYSICO  <- as.character(snakemake@input[["physico"]])[1]

# Outputs
OUT_PDF     <- as.character(snakemake@output[["pdf"]])[1]
OUT_XLSX    <- as.character(snakemake@output[["xlsx"]])[1]
OUT_PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared Plot Parameters (Filtrés pour l'ACP)
PARAM_SHARED            <- snakemake@params[["shared"]]
PARAM_THEME             <- as.character(PARAM_SHARED$theme) %||% "theme_minimal"
PARAM_PDF_SIZE          <- as.numeric(PARAM_SHARED$pdf_size) %||% c(7, 5.5)
PARAM_TITLE_SIZE        <- as.numeric(PARAM_SHARED$title_size) %||% 12
PARAM_SUBTITLE_SIZE     <- as.numeric(PARAM_SHARED$subtitle_size) %||% 10
PARAM_AXES_TITLE_SIZE   <- as.numeric(PARAM_SHARED$axes_title_size) %||% 9.5
PARAM_AXES_TICK_SIZE    <- as.numeric(PARAM_SHARED$axes_tick_size) %||% 8
PARAM_LEGEND_TITLE_SIZE <- as.numeric(PARAM_SHARED$legend_title_size) %||% 8.5
PARAM_LEGEND_SIZE       <- as.numeric(PARAM_SHARED$legend_size) %||% 8
PARAM_PCA_POINT_SIZE    <- as.numeric(PARAM_SHARED$pca_point_size) %||% 2.0
PARAM_PCA_LABEL_SIZE    <- as.numeric(PARAM_SHARED$pca_label_size) %||% 3.0
PARAM_STROKE_WIDTH      <- as.numeric(PARAM_SHARED$stroke_width) %||% 0.5

# Specific PCA Parameters & Translation Templates
TEMPLATE_TITLE    <- as.character(snakemake@params[["title_template"]])[1]
TEMPLATE_SUBTITLE <- as.character(snakemake@params[["subtitle_template"]])[1]
TEXT_SAMPLE_ALL   <- as.character(snakemake@params[["text_sample_all"]])[1] %||% "all samples"

PARAM_RANK        <- tolower(as.character(snakemake@params[["rank"]])[1])
PARAM_TOP_N       <- as.integer(snakemake@params[["top_n"]])[1] %||% 50
PARAM_DIM_X       <- as.integer(snakemake@params[["dim_x"]])[1] %||% 1
PARAM_DIM_Y       <- as.integer(snakemake@params[["dim_y"]])[1] %||% 2
PARAM_PHYSICO_COL <- as.character(snakemake@params[["physico_col"]])
PARAM_STAND_COL   <- as.character(snakemake@params[["stand_col"]])[1]

# Wildcards & Dynamic Variables
WILDCARD_SOURCE   <- tolower(as.character(snakemake@wildcards[["source"]]))[1]
DIM_NAMES         <- paste0("Dim.", c(PARAM_DIM_X, PARAM_DIM_Y))

# ==========================================================================
# Helper : Résolution des textes/titres
# ==========================================================================
resolve_text <- function(template, vars = list()) {
  out <- template
  for (name in names(vars)) {
    out <- gsub(paste0("\\{", name, "\\}"), as.character(vars[[name]]), out)
  }
  return(out)
}

# Dynamic feature_type resolution
FEATURE_TYPE <- if (grepl("level_3", PARAM_RANK, ignore.case = TRUE)) {
  "pathways"
} else if (grepl("gene_description", PARAM_RANK, ignore.case = TRUE) || WILDCARD_SOURCE == "kegg") {
  "genes"
} else {
  "taxa"
}

phyloseq_to_dt_fast <- function(ps) {
  # 1. Extraction matrice otu
  otu_mat <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  
  dt_long <- as.data.table(otu_mat, keep.rownames = "Feature_ID")
  dt_long <- melt(dt_long, id.vars = "Feature_ID", variable.name = "Sample", value.name = PARAM_STAND_COL)
  
  # 2. Métadonnées échantillons
  if (!is.null(phyloseq::sample_data(ps, errorIfNULL = FALSE))) {
    dt_meta <- as.data.table(as(phyloseq::sample_data(ps), "data.frame"), keep.rownames = "Sample")
    dt_long <- merge(dt_long, dt_meta, by = "Sample", all.x = TRUE)
  }
  
  # 3. Taxonomie
  if (!is.null(phyloseq::tax_table(ps, errorIfNULL = FALSE))) {
    dt_tax <- as.data.table(as(phyloseq::tax_table(ps), "matrix"), keep.rownames = "Feature_ID")
    dt_long <- merge(dt_long, dt_tax, by = "Feature_ID", all.x = TRUE)
  }
  
  return(dt_long)
}

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
# Input Files
for (f in c(IN_PHYLOSEQ, IN_PHYSICO)) {
  if (!file.exists(f)) stop(sprintf("❌ Critical Error: Input file does not exist: %s", f))
}

# Rank Parameter
if (is.na(PARAM_RANK) || PARAM_RANK == "" || PARAM_RANK == "null") {
  stop("❌ Critical Error: 'rank' parameter must be configured in config.yaml.")
}

# Dynamic Source Identification
is_kegg   <- grepl("kegg", WILDCARD_SOURCE, ignore.case = TRUE)
is_contig <- grepl("contig", WILDCARD_SOURCE, ignore.case = TRUE)
is_reads  <- grepl("read", WILDCARD_SOURCE, ignore.case = TRUE)

if (!is_kegg && !is_contig && !is_reads) {
  stop(sprintf("❌ Critical Error: Unrecognized source '%s'. Must be kegg, contig, or reads.", WILDCARD_SOURCE))
}

TARGET_ID <- if (is_kegg) "kegg_id" else if (is_contig) "contig_id" else "read_id"

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading Phyloseq object...")
ps <- readRDS(IN_PHYLOSEQ)
if (!inherits(ps, "phyloseq")) stop("❌ Critical Error: Loaded object is not of class 'phyloseq'.")

n_ps_samples <- phyloseq::nsamples(ps)
n_ps_taxa    <- phyloseq::ntaxa(ps)

if (is.null(n_ps_samples) || n_ps_samples == 0) stop("❌ Critical Error: Phyloseq object contains 0 samples.")
if (is.null(n_ps_taxa) || n_ps_taxa == 0) stop("❌ Critical Error: Phyloseq object contains 0 taxa/features.")

message(sprintf("✓ Phyloseq object loaded successfully: %d sample(s) and %d feature(s).", n_ps_samples, n_ps_taxa))

dt_merged <- phyloseq_to_dt_fast(ps)

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PREPARATION
# ------------------------------------------------------------------------------
# Clearing NA entries from the selected RANK
dt_merged[is.na(get(PARAM_RANK)) | get(PARAM_RANK) == "", (PARAM_RANK) := "Unclassified"]

if (is_kegg) {
  feature_col <- "combined_label"

  # Identification dynamique de la colonne KO/kegg_id
  found_col <- intersect("Feature_ID", names(dt_merged))
  
  if (length(found_col) > 0) {
    target_ko_col <- found_col[1]
  } else if ("taxa_id" %in% names(dt_merged)) {
    target_ko_col <- "taxa_id"
  } else {
    target_ko_col <- TARGET_ID
  }

  # Nettoyage et sécurisation des identifiants KO
  ko_vals <- as.character(dt_merged[[target_ko_col]])
  ko_vals[is.na(ko_vals) | ko_vals == "" | ko_vals == "NA"] <- "Unknown_KO"

  # Nettoyage de la description/rang
  rank_vals <- as.character(dt_merged[[PARAM_RANK]])
  rank_vals[is.na(rank_vals) | rank_vals == "" | rank_vals == "Unassigned"] <- "Unknown Function"

  # Concaténation explicite garantissant "KO | Description"
  dt_merged[, (feature_col) := paste(ko_vals, rank_vals, sep = " | ")]

} else {
  feature_col <- PARAM_RANK
}

# Wide-format pivot (Samples x Features)
formula_dcast <- as.formula(paste("date + name ~", feature_col))
tableau_large <- dcast(
  dt_merged,
  formula_dcast,
  value.var = PARAM_STAND_COL,
  fun.aggregate = function(x) sum(x, na.rm = TRUE),
  fill = 0
)

# Safety Replacement of Post-Dcast NAs
n_na <- sum(sapply(tableau_large, function(col) sum(is.na(col))))
if (n_na > 0) message(sprintf("Replacing %d remaining NA(s) with 0 after dcast.", n_na))
for (col in names(tableau_large)) {
  set(tableau_large, which(is.na(tableau_large[[col]])), col, 0)
}

# Loading and Cleaning Physicochemical Data (Excel)
dt_physico <- as.data.table(readxl::read_excel(IN_PHYSICO))
dt_physico[, date := as.character(date)]
dt_physico[, name := as.character(name)]

dt_physico_clean <- dt_physico[, lapply(.SD, function(x) {
  clean_x <- gsub("[^0-9.,-]", "", as.character(x))
  clean_x <- gsub(",", ".", clean_x)
  result  <- as.numeric(clean_x)
  
  n_new_na <- sum(is.na(result) & !is.na(x) & trimws(as.character(x)) != "")
  if (n_new_na > 0) warning(sprintf("%d value(s) could not be parsed as numeric and became NA.", n_new_na))
  
  result
}), .SDcols = PARAM_PHYSICO_COL]

dt_physico_clean[, `:=`(date = dt_physico$date, name = dt_physico$name)]

# Alignment and physicochemical merging with the wide table
tableau_large <- merge(tableau_large, dt_physico_clean, by = c("date", "name"), all.x = TRUE)

# Séparation Métadonnées / Taxonomie / Physico
metadata           <- tableau_large[, .(date = as.factor(date), name = as.factor(name))]
taxon_cols         <- setdiff(names(tableau_large), c("date", "name", PARAM_PHYSICO_COL))
data_pca           <- tableau_large[, taxon_cols, with = FALSE]
physico_for_impute <- as.data.frame(tableau_large[, PARAM_PHYSICO_COL, with = FALSE])

# CLR Transformation & missMDA Allocation
data_pca[is.na(data_pca)] <- 0
donnees_clr <- as.data.frame(as.matrix(compositions::clr(data_pca + 1)))

max_possible_ncp <- min(nrow(physico_for_impute) - 2L, ncol(physico_for_impute) - 1L)
if (max_possible_ncp < 1L) {
  chosen_ncp <- 0L
} else {
  ncp_bound  <- min(2L, max_possible_ncp)
  nb         <- missMDA::estim_ncpPCA(physico_for_impute, ncp.max = ncp_bound)
  chosen_ncp <- nb$ncp
}

res_impute      <- missMDA::imputePCA(physico_for_impute, ncp = chosen_ncp)
physico_complet <- as.data.frame(res_impute$completeObs)

# Running the PCA (FactoMineR)
final_tab_pca  <- cbind(donnees_clr, physico_complet, name = metadata$name)
n_taxons       <- ncol(donnees_clr)
n_physico      <- length(PARAM_PHYSICO_COL)
idx_quanti_sup <- (n_taxons + 1):(n_taxons + n_physico)
idx_quali_sup  <- n_taxons + n_physico + 1L

res_pca <- FactoMineR::PCA(
  final_tab_pca,
  scale.unit = FALSE,
  ncp        = max(PARAM_DIM_X, PARAM_DIM_Y, 5),
  quanti.sup = idx_quanti_sup,
  quali.sup  = idx_quali_sup,
  graph      = FALSE
)

# Extract variance percentage for dynamic axis labeling
var_x <- round(res_pca$eig[PARAM_DIM_X, 2], 1)
var_y <- round(res_pca$eig[PARAM_DIM_Y, 2], 1)

contrib_sum      <- res_pca$var$contrib[, PARAM_DIM_X] + res_pca$var$contrib[, PARAM_DIM_Y]
actual_top_n     <- min(PARAM_TOP_N, length(contrib_sum))
top_taxons       <- names(sort(contrib_sum, decreasing = TRUE)[seq_len(actual_top_n)])
selection_finale <- c(top_taxons, PARAM_PHYSICO_COL)

# ------------------------------------------------------------------------------
# 6. GRAPHICS GENERATION & ANNOTATION
# ------------------------------------------------------------------------------
n_dates   <- nlevels(metadata$date)
n_samples <- nlevels(metadata$name)

color_blind_friendly_base <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
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

# Dynamic titles resolution
title_str <- resolve_text(TEMPLATE_TITLE, list(
  source = toupper(WILDCARD_SOURCE),
  top_n = actual_top_n,
  feature_type = FEATURE_TYPE,
  sample = TEXT_SAMPLE_ALL
))

subtitle_str <- resolve_text(TEMPLATE_SUBTITLE, list(
  dim_x = PARAM_DIM_X,
  dim_y = PARAM_DIM_Y,
  rank = PARAM_RANK
))

# Plot 1: Screeplot Variance
plot_variance <- factoextra::fviz_eig(res_pca, addlabels = TRUE) + 
  get(PARAM_THEME)() +
  theme(
    plot.title     = element_text(size = PARAM_TITLE_SIZE, face = "bold"),
    axis.title     = element_text(size = PARAM_AXES_TITLE_SIZE),
    axis.text      = element_text(size = PARAM_AXES_TICK_SIZE),
    line           = element_line(linewidth = PARAM_STROKE_WIDTH)
  )

# Prepare Individual coordinates data frame for ggplot2 mapping
df_ind_coords <- as.data.frame(res_pca$ind$coord)
df_ind_coords$date <- metadata$date
df_ind_coords$name <- metadata$name

# Plot 2: Biplot PCA
plot_acp <- factoextra::fviz_pca_biplot(
  res_pca,
  axes           = c(PARAM_DIM_X, PARAM_DIM_Y),
  geom.ind       = "none",
  col.var        = "black",
  col.quanti.sup = "blue",
  select.var     = list(name = selection_finale),
  habillage      = "none",
  mean.point     = FALSE,
  repel          = TRUE,
  repel.opts     = list(max.overlaps = Inf, force = 2),
  labelsize      = PARAM_PCA_LABEL_SIZE,
  title          = title_str,
  subtitle       = subtitle_str
) +
  geom_point(
    data = df_ind_coords,
    aes(
      x     = .data[[DIM_NAMES[1]]],
      y     = .data[[DIM_NAMES[2]]],
      shape = date,
      color = name
    ),
    size = PARAM_PCA_POINT_SIZE, alpha = 0.8
  ) +
  scale_shape_manual(
    values = shapes_dates[seq_len(n_dates)],
    name   = "Dates",
  ) +
  scale_color_manual(
    values = palette_samples[seq_len(n_samples)],
    name   = "Samples",
  ) +
  labs(
    x = sprintf("Dim %d (%s%%)", PARAM_DIM_X, var_x),
    y = sprintf("Dim %d (%s%%)", PARAM_DIM_Y, var_y)
  ) +
  get(PARAM_THEME)() +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = PARAM_TITLE_SIZE, face = "bold"),
    plot.subtitle   = element_text(size = PARAM_SUBTITLE_SIZE),
    axis.title      = element_text(size = PARAM_AXES_TITLE_SIZE),
    axis.text       = element_text(size = PARAM_AXES_TICK_SIZE),
    legend.title    = element_text(size = PARAM_LEGEND_TITLE_SIZE, face = "bold"),
    legend.text     = element_text(size = PARAM_LEGEND_SIZE),
    line            = element_line(linewidth = PARAM_STROKE_WIDTH)
  )

# Contribution Table
contrib_top <- as.data.table(res_pca$var$contrib[top_taxons, , drop = FALSE], keep.rownames = "Taxon")
table_grob  <- gridExtra::tableGrob(contrib_top)

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
message("INFO: Generating PCA graphics PDF...")

with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  pages_rendered <- 0
  
  if (!is.null(plot_variance)) {
    render_page(plot_variance)
    pages_rendered <- pages_rendered + 1
  }
  
  if (!is.null(plot_acp)) {
    render_page(plot_acp)
    pages_rendered <- pages_rendered + 1
  }
  
  if (!is.null(table_grob)) {
    # Render table smoothly inside standard ggplot device context
    render_page(ggplot2::qplot(1:10, 1:10, geom = "blank") + 
                  annotation_custom(table_grob) + 
                  theme_void())
    pages_rendered <- pages_rendered + 1
  }
  
  if (pages_rendered == 0) {
    render_fallback("⚠️ WARNING: Empty PCA output: No plots were generated.")
  }
})

# Export Parquet Contact Information for Individuals
coords_ind <- as.data.table(res_pca$ind$coord)
coords_ind[, name := as.character(metadata$name)]
coords_ind[, date := as.character(metadata$date)]
coords_ind[, (PARAM_PHYSICO_COL) := physico_complet]

if (nrow(coords_ind) > 0) {
  writexl::write_xlsx(contrib_top, OUT_XLSX)
  arrow::write_parquet(coords_ind, OUT_PARQUET)

  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Excel   : ", OUT_XLSX)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  writexl::write_xlsx(data.frame(), OUT_XLSX)
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: Empty output: No results written.")
}