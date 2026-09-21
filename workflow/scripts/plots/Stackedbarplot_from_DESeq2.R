# ==============================================================================
# PROJECT : MicrobExplorer
# SCRIPT  : Stackedbarplot from DESeq2 Results
# PURPOSE : Generate stacked barplots for DESeq2 differentially abundant features
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
  library(viridis)
  library(arrow)
  library(rlang)
  library(glue)
  # Libraries Bioconductor
  library(phyloseq)
})

# Load helper
helper_path <- dirname(snakemake@scriptdir)
source(file.path(helper_path, "utils", "utils_stackedbarplot_from_DESeq2.R"))
source("workflow/scripts/utils/utils_pdf.R")

# Disable automatic factors and set strict mode
options(stringsAsFactors = FALSE, warn = 1)

# ------------------------------------------------------------------------------
# 2. SNAKEMAKE I/O & PARAMETERS BINDING
# ------------------------------------------------------------------------------
# Inputs
IN_DESEQ    <- as.character(snakemake@input[["deseq_files"]])[1]
IN_PHYLOSEQ <- as.character(snakemake@input[["phyloseq_obj"]])[1]

# Outputs
OUT_PDF     <- as.character(snakemake@output[["pdf"]])[1]
OUT_PARQUET <- as.character(snakemake@output[["parquet"]])[1]

# Shared Plot Parameters
PARAM_SHARED   <- snakemake@params[["shared"]]
PARAM_PDF_SIZE <- c(12, 10)  # as.numeric(PARAM_SHARED$pdf_size) %||% c(7, 5.5)
PARAM_PALETTE  <- tolower(as.character(PARAM_SHARED$palette)) %||% "turbo"

# Specific Stackedbarplot DESeq2 Parameters & Translation Templates
TEMPLATE_TITLE    <- as.character(snakemake@params[["title_template"]])[1]
TEMPLATE_SUBTITLE <- as.character(snakemake@params[["subtitle_template"]])[1]
TEXT_SAMPLE_ALL   <- as.character(snakemake@params[["text_sample_all"]])[1] %||% "all samples"

PARAM_CONTRAST <- unlist(as.character(snakemake@params[["contrast"]]))
PARAM_PADJ     <- as.numeric(unlist(snakemake@params[["padj"]])) %||% c(0.05, 0.05, 0.05)
PARAM_LFC      <- as.numeric(unlist(snakemake@params[["lfc"]]))  %||% c(1.3, 1.3, 1.3)
PARAM_TOP_N    <- as.integer(snakemake@params[["top_n"]])[1] %||% 10
PARAM_RANK     <- tolower(as.character(snakemake@params[["rank"]]))[1]

# Wildcards & Dynamic Variables
WILDCARD_SOURCE <- tolower(as.character(snakemake@wildcards[["source"]]))[1]

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

# ------------------------------------------------------------------------------
# 3. PARAMETER VALIDATION ("FAIL-FAST")
# ------------------------------------------------------------------------------
# Validate the required files
for (f in c(IN_DESEQ, IN_PHYLOSEQ)) {
  if (!file.exists(f)) stop(sprintf("❌ Critical Error: Input file does not exist: %s", f))
}

# Confirm Rank
if (is.null(PARAM_RANK) || is.na(PARAM_RANK) || PARAM_RANK == "" || PARAM_RANK == "null") {
  stop("❌ Critical Error: 'rank' parameter must be configured in config.yaml.")
}

# Identify the source
is_kegg   <- grepl("kegg", WILDCARD_SOURCE, ignore.case = TRUE)
is_contig <- grepl("contig", WILDCARD_SOURCE, ignore.case = TRUE)
is_reads  <- grepl("read", WILDCARD_SOURCE, ignore.case = TRUE)

if (!is_reads && !is_contig && !is_kegg) {
  stop(sprintf("❌ Critical Error: Unrecognized source '%s'. Must be reads, contig, or kegg", WILDCARD_SOURCE))
}

TARGET_ID <- if (is_kegg) "kegg_id" else if (is_contig) "contig_id" else "read_id"

valid_contrasts <- c("ref", "combo", "date")

if (length(PARAM_CONTRAST) == 0 || 
    any(is.na(PARAM_CONTRAST)) || 
    !all(PARAM_CONTRAST %in% valid_contrasts)) {
  
  stop(sprintf(
    "❌ Critical Error: Invalid contrast parameter(s) [%s]. All values must be among: %s",
    paste(PARAM_CONTRAST, collapse = ", "),
    paste(valid_contrasts, collapse = ", ")
  ))
}

# ------------------------------------------------------------------------------
# 4. DATA LOADING & INTEGRITY CHECKS
# ------------------------------------------------------------------------------
message("INFO: Loading Phyloseq object...")
ps <- readRDS(IN_PHYLOSEQ)

if (!inherits(ps, "phyloseq")) stop("❌ Critical Error: Loaded object is not of class 'phyloseq'.")
if (phyloseq::nsamples(ps) == 0 || phyloseq::ntaxa(ps) == 0) stop("❌ Phyloseq object is empty.")

message("INFO: Loading DESeq2 Master table...")
master_rds <- readRDS(IN_DESEQ)
master_table <- if (is.list(master_rds) && !is.null(master_rds$dt)) {
  rbindlist(master_rds$dt, use.names = TRUE, fill = TRUE)
} else {
  as.data.table(master_rds)
}

if (!"Feature_ID" %in% names(master_table)) {
  if ("KO_Number" %in% names(master_table)) setnames(master_table, "KO_Number", "Feature_ID")
  else setnames(master_table, 1, "Feature_ID")
}

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATIONS & PREPARATION
# ------------------------------------------------------------------------------
# Relativization and Conversion to data.table
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
df_abundance <- phyloseq_to_dt(ps_rel)

if (!"OTU" %in% names(df_abundance)) {
  if (TARGET_ID %in% names(df_abundance)) {
    setnames(df_abundance, TARGET_ID, "OTU")
  } else {
    setnames(df_abundance, 1, "OTU")
  }
}
df_abundance[, (TARGET_ID) := OTU]

if ("gene_description" %in% names(master_table) || "Description" %in% names(master_table)) {
  desc_col <- if ("gene_description" %in% names(master_table)) "gene_description" else "Description"
  
  # Extraire une table de correspondance unique Feature_ID -> Description
  gene_map <- unique(master_table[!is.na(Feature_ID), .(OTU = Feature_ID, gene_description = get(desc_col))])
  gene_map <- gene_map[!duplicated(OTU)]
  
  # Nettoyage des numéros EC : retire [EC:x.x.x.x] ou (EC:x.x.x.x)
  gene_map[, gene_description := trimws(gsub("\\[EC:.*?\\]|\\(EC:.*?\\)", "", gene_description))]
  
  # Fusionner avec df_abundance
  df_abundance <- merge(df_abundance, gene_map, by = "OTU", all.x = TRUE)
}

# Helper pour formater la description des gènes (Premier gène + (+ n ...))
format_gene_desc <- function(raw_val) {
  if (is.na(raw_val) || raw_val == "" || raw_val == "NA" || raw_val == "Unassigned") return(NA_character_)
  
  # 1. Découpage des blocs gènes séparés par une VIRGULE
  gene_blocks <- unlist(strsplit(as.character(raw_val), ","))
  
  # 2. Nettoyage de chaque bloc : suppression de ce qui suit le POINT-VIRGULE
  clean_genes <- sapply(gene_blocks, function(block) {
    trimws(sub(";.*$", "", block))
  })
  
  # Conservation des noms non vides
  clean_genes <- clean_genes[clean_genes != ""]
  
  # 3. Assemblage du label et décompte des gènes restants
  if (length(clean_genes) > 1) {
    return(sprintf("%s (+ %d ...)", clean_genes[1], length(clean_genes) - 1))
  } else if (length(clean_genes) == 1) {
    return(clean_genes[1])
  } else {
    return(raw_val)
  }
}

# Abundance filtering based on the tested features
all_tested_features <- unique(master_table$Feature_ID)
df_abundance <- df_abundance[OTU %in% all_tested_features]
setkey(df_abundance, OTU)

results_list  <- split(master_table, master_table$Comparison)
all_plot_data <- list()

# ------------------------------------------------------------------------------
# 6. PROCESSING & VISUALIZATION LOOP
# ------------------------------------------------------------------------------
with_pdf(OUT_PDF, PARAM_PDF_SIZE, {
  page_count <- 0

  for (i in seq_along(PARAM_CONTRAST)) {
    contrast_type <- PARAM_CONTRAST[i]
    current_padj  <- PARAM_PADJ[i]
    current_lfc   <- PARAM_LFC[i]
    
    comparisons_this_type <- unique(master_table[Contrast_Type == contrast_type]$Comparison)
    if (length(comparisons_this_type) == 0) {
      clean_contrast <- gsub("^group_", "", contrast_type)
      comparisons_this_type <- unique(master_table[gsub("^group_", "", Comparison) == clean_contrast]$Comparison)
    }
    
    if (length(comparisons_this_type) == 0) next
    
    for (comp in comparisons_this_type) {
      deseq_result <- results_list[[comp]]
      if (is.null(deseq_result)) next
      
      # Process Direction: UP / DOWN
      for (dir in c("up", "down")) {
        sig_features <- get_significant_features(deseq_result, current_padj, current_lfc, dir)
        if (length(sig_features) == 0) next
        
        df_plot <- prepare_plot_data(
          df_abundance, 
          sig_features, 
          TOP_N = PARAM_TOP_N,
          group_column = "group", 
          rank_column = PARAM_RANK,
          comp_name = comp, 
          time_column = "date"
        )
        
        cats <- setdiff(levels(df_plot$Category), c("Not significant", "Significant (Other)"))
        colors <- generate_color_palette(cats, palette_name = PARAM_PALETTE)
        
        actual_top_n <- min(PARAM_TOP_N, length(cats))

        title_str <- resolve_text(TEMPLATE_TITLE, list(
          source = toupper(WILDCARD_SOURCE),
          top_n = actual_top_n,
          feature_type = FEATURE_TYPE,
          contrast = toupper(comp)
        ))

        subtitle_str <- resolve_text(TEMPLATE_SUBTITLE, list(
          padj_threshold = current_padj,
          lfc_threshold = current_lfc,
          rank = PARAM_RANK
        ))

        resolved_title <- paste0(title_str, " | ", toupper(dir), "-REGULATED")
        
        p <- create_stackedbarplot(df_plot, title = resolved_title, subtitle = subtitle_str, feature_colors = colors, group_label = "group")
        
        render_page(p)
        page_count <- page_count + 1

        df_plot[, `:=`(Comparison = comp, Contrast_Type = contrast_type, Direction = dir)]
        all_plot_data[[paste0(comp, "_", dir)]] <- df_plot

      }
    }
  }

  if (page_count == 0) {
    render_fallback("No significant features found for Stackedbar plots.")
  }
})

# ------------------------------------------------------------------------------
# 7. EXPORTS & OUTPUT GENERATION
# ------------------------------------------------------------------------------
if (length(all_plot_data) > 0) {
  final_dt <- rbindlist(all_plot_data, fill = TRUE)
  arrow::write_parquet(final_dt, OUT_PARQUET)
  message("✓ Success exports written:")
  message("  - PDF     : ", OUT_PDF)
  message("  - Parquet : ", OUT_PARQUET)
} else {
  arrow::write_parquet(data.table(), OUT_PARQUET)
  warning("⚠️ WARNING: Empty output: No significant results found.")
}