################################################################################
# Project : "MicrobExplorer"
# Script: "Stackedbarplot from DESeq2 Results"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# Libraries CRAN
library(data.table)
library(ggplot2)
library(viridis)
library(arrow)
library(rlang)
library(glue)

# Libraries Bioconductor
library(phyloseq)

# Load helper
helper_path <- dirname(snakemake@scriptdir)
source(file.path(helper_path, "utils", "utils_stackedbarplot_from_DESeq2.R"))

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================

# Inputs
DESEQ_FILES <- as.character(snakemake@input[["deseq_files"]])[1]
PHYLOSEQ_OBJ <- as.character(snakemake@input[["phyloseq_obj"]])[1]
METADATA <- as.character(snakemake@input[["metadata"]])[1]

# Outputs
PDF <- as.character(snakemake@output[["pdf"]])[1]
PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared plots features
SHARED <- snakemake@params[["shared"]]
PDF_SIZE    <- as.numeric(SHARED$pdf_size) %||% c(14, 12)
PALETTE   <- tolower(as.character(SHARED$palette)) %||% "turbo"

# Parameters
TITLE_TEMPLATE    <- as.character(snakemake@params[["title"]])[1] %||% "{source} | Top {top_n} {feature_type} differentially abundant in {contrast}"
SUBTITLE_TEMPLATE <- as.character(snakemake@params[["subtitle"]])[1] %||% "Pvalue threshold : {padj_threshold} | Log2 Fold Change threshold : {lfc_threshold} | {rank}"
CONTRAST_LIST <- unlist(as.character(snakemake@params[["contrast"]]))
PADJ_THRESHOLD <- as.numeric(unlist(snakemake@params[["padj"]])) %||% c(0.05, 0.05, 0.05)
LFC_THRESHOLD  <- as.numeric(unlist(snakemake@params[["lfc"]]))  %||% c(1.3, 1.3, 1.3)
TOP_N <- as.integer(snakemake@params[["top_n"]])[1] %||% 10
RANK <- tolower(as.character(snakemake@params[["rank"]]))[1]

# Wildcards
SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

# Standardization of Structure Type
FEATURE_TYPE <- if (SOURCE == "kegg") "pathways" else "taxons"
# ==========================================================================
# 1. Data Loading
# ==========================================================================
message("Loading precomputed data from Master RDS...")
# 🟢 ADAPTATION : Lecture du slot $dt depuis le nouveau Master RDS
master_rds_obj <- readRDS(DESEQ_FILES)

if (is.list(master_rds_obj) && !is.null(master_rds_obj$dt)) {
  # Si master_rds_obj$dt est une liste (ref, date, combo), on combine tout en une seule master_table
  master_table <- rbindlist(master_rds_obj$dt, use.names = TRUE, fill = TRUE)
} else {
  master_table <- as.data.table(master_rds_obj)
}

# 🟢 ADAPTATION : Standardisation du nom de la colonne d'identifiant
if (!"Feature_ID" %in% names(master_table) && "KO_Number" %in% names(master_table)) {
  setnames(master_table, "KO_Number", "Feature_ID")
} else if (!"Feature_ID" %in% names(master_table)) {
  setnames(master_table, 1, "Feature_ID")
}

results_list <- split(master_table, master_table$Comparison)
message(sprintf("  → Master table compiled: %d comparisons across %s contrast types", 
                length(results_list), paste(unique(master_table$Contrast_Type), collapse = ", ")))

# Load and prepare phyloseq object
ps <- readRDS(PHYLOSEQ_OBJ)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
df_abundance <- phyloseq_to_dt(ps_rel)
message(sprintf("✓ Loaded phyloseq object: %d OTUs, %d samples", nrow(otu_table(ps)), ncol(otu_table(ps))))

# ✅ TON OPTIMISATION PRÉSERVÉE : réduire df_abundance aux seules features testées par DESeq2
all_tested_features <- unique(master_table$Feature_ID)
n_before <- uniqueN(df_abundance$OTU)
df_abundance <- df_abundance[OTU %in% all_tested_features]
message(sprintf("  → df_abundance réduit de %d à %d OTU (features testées par DESeq2 uniquement)", 
                n_before, uniqueN(df_abundance$OTU)))

# ✅ TON OPTIMISATION PRÉSERVÉE : indexer sur OTU pour accélérer les filtres %in% répétés
setkey(df_abundance, OTU)

# Accumulateur pour l'export parquet final
all_plot_data <- list()

pdf(PDF, width = PDF_SIZE[1], height = PDF_SIZE[2])

for (i in seq_along(CONTRAST_LIST)) {
  contrast_type <- CONTRAST_LIST[i]
  CURRENT_PADJ  <- PADJ_THRESHOLD[i]
  CURRENT_LFC   <- LFC_THRESHOLD[i]
  
  # Récupération des comparaisons associées à ce type de contraste
  comparisons_this_type <- unique(master_table[Contrast_Type == contrast_type]$Comparison)
  
  # Fallback si l'utilisateur a écrit "group_TD2_vs_TD1" au lieu de "TD2_vs_TD1" ou inversement
  if (length(comparisons_this_type) == 0) {
    clean_contrast <- gsub("^group_", "", contrast_type)
    comparisons_this_type <- unique(master_table[gsub("^group_", "", Comparison) == clean_contrast]$Comparison)
  }
  
  if (length(comparisons_this_type) == 0) {
    message(sprintf("[WARNING] Aucun Comparison trouvé pour Contrast_Type '%s' (position %d). Skip.", contrast_type, i))
    next
  }
  
  message(sprintf("=== Contrast_Type '%s' : %d comparisons, padj <= %g, lfc >= %g ===", 
                  contrast_type, length(comparisons_this_type), CURRENT_PADJ, CURRENT_LFC))
  
  for (comp in comparisons_this_type) {
    deseq_result <- results_list[[comp]]
    if (is.null(deseq_result)) next
    
    message(sprintf("  Processing comparison: %s", comp))
    
    # --- A. Up-regulated ---
    sig_up <- get_significant_features(deseq_result, CURRENT_PADJ, CURRENT_LFC, "up")
    if (length(sig_up) > 0) {
      message(sprintf("    → %d up-regulated features", length(sig_up)))
      df_plot_up <- prepare_plot_data(
        df_abundance, 
        sig_up, 
        TOP_N = TOP_N,
        group_column = "group", 
        rank_column = RANK,
        comp_name    = comp,
        time_column  = "date"
      )
      
      current_top_cats <- levels(df_plot_up$Category)
      current_top_cats <- current_top_cats[!current_top_cats %in% c("Not significant", "Significant (Other)")]
      my_feature_colors <- generate_color_palette(current_top_cats, palette_name = PALETTE)
      
      resolved_title_up <- glue(TITLE_TEMPLATE, source=toupper(SOURCE), top_n = TOP_N, feature_type = FEATURE_TYPE, contrast = toupper(comp), " | UP-REGULATED")
      resolved_subtitle_up <- glue(SUBTITLE_TEMPLATE, padj_threshold = CURRENT_PADJ, lfc_threshold = CURRENT_LFC, rank = RANK)
      
      p_up <- create_stackedbarplot(
        df_plot_up, title = resolved_title_up, subtitle = resolved_subtitle_up,
        feature_colors = my_feature_colors, group_label = "group"
      )
      print(p_up)
      
      df_plot_up[, `:=`(Comparison = comp, Contrast_Type = contrast_type, Direction = "up")]
      all_plot_data[[paste0(comp, "_up")]] <- df_plot_up
    }
    
    # --- B. Down-regulated ---
    sig_down <- get_significant_features(deseq_result, CURRENT_PADJ, CURRENT_LFC, "down")
    if (length(sig_down) > 0) {
      message(sprintf("    → %d down-regulated features", length(sig_down)))
      df_plot_down <- prepare_plot_data(
        df_abundance, sig_down, 
        TOP_N = TOP_N,
        group_column = "group", 
        rank_column = RANK,
        comp_name    = comp,          # <--- AJOUT : La variable de ta boucle (ex: "J15_vs_J1")
        time_column  = "date"
      )
      
      current_top_cats <- levels(df_plot_down$Category)
      current_top_cats <- current_top_cats[!current_top_cats %in% c("Not significant", "Significant (Other)")]
      my_feature_colors <- generate_color_palette(current_top_cats, palette_name = PALETTE)
      
      resolved_title_down <- glue(TITLE_TEMPLATE, source=toupper(SOURCE), top_n = TOP_N, feature_type = FEATURE_TYPE, contrast = toupper(comp), " | DOWN-REGULATED")
      resolved_subtitle_down <- glue(SUBTITLE_TEMPLATE, padj_threshold = CURRENT_PADJ, lfc_threshold = CURRENT_LFC, rank = RANK)
      
      p_down <- create_stackedbarplot(
        df_plot_down, title = resolved_title_down, subtitle = resolved_subtitle_down,
        feature_colors = my_feature_colors, group_label = "group"
      )
      print(p_down)
      
      df_plot_down[, `:=`(Comparison = comp, Contrast_Type = contrast_type, Direction = "down")]
      all_plot_data[[paste0(comp, "_down")]] <- df_plot_down
    }
    
    if (length(sig_up) == 0 && length(sig_down) == 0) {
      message(sprintf("    → No significant features found with padj < %g and lfc > %g", CURRENT_PADJ, CURRENT_LFC))
    }
  }
}
dev.off()

# Export final du parquet attendu par Snakemake
if (length(all_plot_data) > 0) {
  final_dt <- rbindlist(all_plot_data, fill = TRUE)
  export_results(final_dt, PARQUET)
} else {
  warning("Aucune donnée significative trouvée : le parquet de sortie sera vide ou non créé.")
}